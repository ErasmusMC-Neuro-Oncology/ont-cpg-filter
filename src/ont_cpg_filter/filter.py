"""Core logic: keep only 5mC/5hmC calls on CpG positions in MM/ML tags.

MM format (one entry per canonical base and strand), e.g. ``C+mh?,0,3;``:
  - base + strand + modification codes + optional flag ('?' or '.')
  - one skip count per call: the number of unmodified candidate bases
    skipped before that call
ML holds one probability per call and code, in the same order.

Orientation: MM refers to the original read. For reverse-aligned reads the
stored SEQ is its reverse complement, so the strand flips and the candidate
order runs from the end of SEQ.

The CpG decision is made only from the reference (see Reference.site_mask).
SEQ is used only to count candidate bases for the skip counts.
"""

from __future__ import annotations

import logging
import time
from array import array
from dataclasses import dataclass, field

import numpy as np
import pysam

from ont_cpg_filter.reference import Reference

log = logging.getLogger(__name__)

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
PROGRESS_EVERY = 1_000_000


@dataclass
class ModEntry:
    base: str  # canonical base as written in MM: 'C' or 'G'
    strand: str  # '+' or '-'
    codes: str  # modification codes, e.g. 'm' or 'mh'
    flag: str  # '', '?' or '.'
    skips: list[int] = field(default_factory=list)

    def to_mm(self) -> str:
        skips = "".join(f",{s}" for s in self.skips)
        return f"{self.base}{self.strand}{self.codes}{self.flag}{skips};"


def parse_mm(mm: str) -> list[ModEntry]:
    entries = []
    for part in mm.split(";"):
        if not part:
            continue
        fields = part.split(",")
        header = fields[0]
        if len(header) < 3 or header[1] not in "+-":
            raise ValueError(f"malformed MM entry: {part!r}")
        rest = header[2:]
        flag = ""
        if rest and rest[-1] in "?.":
            flag, rest = rest[-1], rest[:-1]
        if not rest.isalpha():
            raise NotImplementedError(
                f"numeric ChEBI modification codes are not supported: {part!r}"
            )
        skips = [int(x) for x in fields[1:] if x != ""]
        entries.append(ModEntry(header[0], header[1], rest, flag, skips))
    return entries


def serialize_mm(entries: list[ModEntry]) -> str:
    return "".join(e.to_mm() for e in entries)


@dataclass
class Stats:
    reads: int = 0
    reads_with_mm: int = 0
    kept: int = 0
    removed: int = 0


def filter_read(read: pysam.AlignedSegment, ref: Reference, stats: Stats) -> None:
    """Remove modification calls whose reference position is not a CpG site.

    Calls that are removed make their entry's flag '?' (unknown), not '.'
    (unmodified), so the removed calls are not reported as unmodified.
    Calls on soft-clipped bases (no reference position) are removed too.
    """
    seq = read.query_sequence
    if seq is None:
        return

    entries = parse_mm(read.get_tag("MM"))
    ml_in = list(read.get_tag("ML")) if read.has_tag("ML") else None

    chrom = read.reference_name
    seq_arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    # Reference position per SEQ index; -1 for soft clips.
    rpos = np.fromiter(
        (-1 if r is None else r for r in read.get_reference_positions(full_length=True)),
        dtype=np.int64,
        count=len(seq),
    )

    new_entries: list[ModEntry] = []
    ml_out: list[int] = []
    ml_pos = 0

    for entry in entries:
        # MM refers to the original read. For reverse-aligned reads SEQ is its
        # reverse complement, so the strand flips relative to the stored SEQ.
        strand = entry.strand
        if read.is_reverse:
            strand = "-" if strand == "+" else "+"
        seq_base = entry.base if strand == "+" else entry.base.translate(_COMPLEMENT)

        # SEQ indices of candidate bases, in the order the skip counts run.
        cand = np.flatnonzero(seq_arr == ord(seq_base))
        if strand == "-":
            cand = cand[::-1]
        n_codes = len(entry.codes)

        # Skip counts -> indices into cand: call k = cumsum(skip + 1) - 1.
        skips = np.asarray(entry.skips, dtype=np.int64)
        calls = np.cumsum(skips + 1) - 1
        n_calls = calls.size
        if n_calls and calls[-1] >= cand.size:
            raise ValueError(
                f"{read.query_name}: MM entry {entry.base}{entry.strand}{entry.codes} needs "
                f"{calls[-1] + 1} {seq_base} bases but SEQ has {cand.size} "
                f"(seq length {len(seq)}, MM={read.get_tag('MM')!r})"
            )

        n_vals = n_calls * n_codes
        vals = ml_in[ml_pos : ml_pos + n_vals] if ml_in is not None else None
        ml_pos += n_vals

        if n_calls:
            keep = ref.site_mask(
                chrom,
                rpos[cand[calls]],
                forward=(seq_base == ref.motif[0]),
            )
            kept = np.flatnonzero(keep)
        else:
            kept = np.zeros(0, dtype=np.int64)

        removed = n_calls - kept.size
        stats.kept += kept.size * n_codes
        stats.removed += removed * n_codes

        if kept.size == 0:
            continue

        kept_k = calls[kept]
        entry.skips = (np.diff(np.concatenate(([-1], kept_k))) - 1).tolist()
        if removed:
            entry.flag = "?"
        new_entries.append(entry)

        if vals is not None:
            for j in kept.tolist():
                ml_out.extend(vals[j * n_codes : (j + 1) * n_codes])

    if new_entries:
        read.set_tag("MM", serialize_mm(new_entries))
        if ml_in is not None:
            read.set_tag("ML", array("B", ml_out))
    else:
        read.set_tag("MM", None)
        if read.has_tag("ML"):
            read.set_tag("ML", None)


def filter_bam(
    in_path: str,
    out_path: str,
    fasta_path: str,
    motif: str = "CG",
    threads: int = 1,
) -> Stats:
    stats = Stats()
    start = time.monotonic()
    log.info("input: %s", in_path)
    log.info("output: %s", out_path)
    log.info("motif: %s", motif)
    with pysam.AlignmentFile(in_path, "rb", threads=threads) as inbam, \
            Reference(fasta_path, motif=motif) as ref, \
            pysam.AlignmentFile(out_path, "wb", template=inbam, threads=threads) as outbam:
        log.info("processing reads")
        for read in inbam.fetch(until_eof=True):
            stats.reads += 1
            if not read.is_unmapped and read.has_tag("MM"):
                stats.reads_with_mm += 1
                filter_read(read, ref, stats)
            outbam.write(read)
            if stats.reads % PROGRESS_EVERY == 0:
                log.info(
                    "%d reads processed (%d with MM), %.0f s elapsed",
                    stats.reads, stats.reads_with_mm, time.monotonic() - start,
                )
    log.info(
        "done: %d reads (%d with MM) in %.1f s",
        stats.reads, stats.reads_with_mm, time.monotonic() - start,
    )
    return stats
