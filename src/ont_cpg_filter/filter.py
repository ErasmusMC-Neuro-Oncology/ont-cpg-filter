"""Core logic: keep only 5mC/5hmC calls on CpG positions in MM/ML tags.

MM format (one entry per canonical base and strand), e.g. ``C+mh?,0,3;``:
  - base + strand + modification codes + optional flag ('?' or '.')
  - one skip count per call: the number of candidate bases skipped before
    that call (candidates are all C's, or all G's for the '-' strand, in SEQ)
ML holds one probability per call and code, in the same order.

Output rule:
  - every CpG position in the reference gets an explicit call. Its ML value is
    kept if the input had a call there, otherwise 0 (unmodified). This keeps
    the meaning the input had for unlisted CpG C's under flag '.'.
  - every other position is unlisted, and the flag is '?' (no information),
    so it is not coloured in IGV.

Orientation: MM refers to the original read. For reverse-aligned reads the
stored SEQ is its reverse complement, so the strand flips and the candidate
order runs from the end of SEQ.

The CpG decision is made only from the reference (Reference.site_mask).
SEQ is used only to find candidate bases.
"""

from __future__ import annotations

import logging
import os
import time
from array import array
from dataclasses import dataclass, field

import numpy as np
import pysam
from tqdm import tqdm

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
    skips: np.ndarray | list = field(default_factory=list)

    def to_mm(self) -> str:
        skips = ",".join(map(str, self.skips))
        sep = "," if len(self.skips) else ""
        return f"{self.base}{self.strand}{self.codes}{self.flag}{sep}{skips};"


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
        # Parse all skip counts at once in C instead of one int() per element.
        skips = np.fromstring(",".join(fields[1:]), dtype=np.int64, sep=",")
        entries.append(ModEntry(header[0], header[1], rest, flag, skips))
    return entries


def serialize_mm(entries: list[ModEntry]) -> str:
    return "".join(e.to_mm() for e in entries)


@dataclass
class Stats:
    reads: int = 0
    reads_with_mm: int = 0
    kept: int = 0  # existing calls on CpG positions that were kept
    removed: int = 0  # existing calls outside CpG that were removed
    added: int = 0  # CpG positions without a call, added with ML=0


def count_reads_from_index(bam_path: str) -> int | None:
    """Total number of reads from the BAM index (no pass over the file).

    Returns None if the BAM has no .bai or .csi index.
    """
    if not (os.path.exists(bam_path + ".bai") or os.path.exists(bam_path + ".csi")):
        return None
    total = 0
    for line in pysam.idxstats(bam_path).splitlines():
        parts = line.split("\t")
        if len(parts) >= 4:  # ref, length, mapped, unmapped
            total += int(parts[2]) + int(parts[3])
    return total


def ref_blocks(read: pysam.AlignedSegment) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(query start, reference start, length) per M/=/X CIGAR block, sorted by query start."""
    q_starts, r_starts, lengths = [], [], []
    qpos = 0
    rpos = read.reference_start
    for op, n in read.cigartuples or ():
        if op in (0, 7, 8):  # M, =, X
            q_starts.append(qpos)
            r_starts.append(rpos)
            lengths.append(n)
            qpos += n
            rpos += n
        elif op in (1, 4):  # I, S
            qpos += n
        elif op in (2, 3):  # D, N
            rpos += n
    return (
        np.asarray(q_starts, dtype=np.int64),
        np.asarray(r_starts, dtype=np.int64),
        np.asarray(lengths, dtype=np.int64),
    )


def query_to_ref_at(blocks: tuple[np.ndarray, np.ndarray, np.ndarray], qidx: np.ndarray) -> np.ndarray:
    """Reference position for the given SEQ indices (-1 for soft clips and insertions).

    Binary search over the CIGAR blocks: the cost depends on how many indices
    are queried, not on the read length.
    """
    q0, r0, ln = blocks
    out = np.full(qidx.shape, -1, dtype=np.int64)
    if q0.size == 0 or qidx.size == 0:
        return out
    b = np.searchsorted(q0, qidx, side="right") - 1
    bb = np.maximum(b, 0)
    inside = (b >= 0) & (qidx < q0[bb] + ln[bb])
    out[inside] = (r0[bb] + (qidx - q0[bb]))[inside]
    return out


def query_to_ref(read: pysam.AlignedSegment, seq_len: int) -> np.ndarray:
    """Reference position per SEQ index (-1 for soft clips and insertions). Full array."""
    return query_to_ref_at(ref_blocks(read), np.arange(seq_len, dtype=np.int64))


def filter_read(read: pysam.AlignedSegment, ref: Reference, stats: Stats) -> None:
    """Rewrite MM/ML so that only CpG positions carry a call.

    CpG positions: existing call kept (its ML value), or added with ML=0.
    Other positions: not listed, flag '?'.
    Soft-clipped and inserted bases have no reference position and are never CpG.
    """
    seq = read.query_sequence
    if seq is None:
        return

    entries = parse_mm(read.get_tag("MM"))
    ml_in = list(read.get_tag("ML")) if read.has_tag("ML") else None

    chrom = read.reference_name
    seq_arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    blocks = ref_blocks(read)

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

        # SEQ indices of all candidate bases, in the order the skip counts run.
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

        # Which candidate bases are CpG sites, decided by the reference only.
        cand_site = ref.site_mask(
            chrom,
            query_to_ref_at(blocks, cand),
            forward=(seq_base == ref.motif[0]),
        )
        kept_j = np.flatnonzero(cand_site[calls]) if n_calls else np.zeros(0, dtype=np.int64)
        row_of = {int(calls[j]): int(j) for j in kept_j}  # cand index -> row in ML
        cpg_idx = np.flatnonzero(cand_site)

        final = np.union1d(calls[kept_j], cpg_idx).astype(np.int64)

        stats.kept += kept_j.size * n_codes
        stats.removed += (n_calls - kept_j.size) * n_codes
        stats.added += (final.size - kept_j.size) * n_codes

        if final.size == 0:
            continue

        entry.skips = (np.diff(np.concatenate(([-1], final))) - 1).tolist()
        entry.flag = "?"
        new_entries.append(entry)

        if ml_in is not None:
            for k in final.tolist():
                j = row_of.get(k)
                if j is not None:
                    ml_out.extend(vals[j * n_codes : (j + 1) * n_codes])
                else:
                    ml_out.extend([0] * n_codes)

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
    index_output: bool = True,
    show_progress: bool = False,
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
        # Progress bar only when asked for: update() per read costs time.
        bar = (
            tqdm(total=count_reads_from_index(in_path), unit="read")
            if show_progress
            else None
        )
        for read in inbam.fetch(until_eof=True):
            stats.reads += 1
            if bar is not None:
                bar.update(1)
            if not read.is_unmapped and read.has_tag("MM"):
                stats.reads_with_mm += 1
                filter_read(read, ref, stats)
            outbam.write(read)
            if stats.reads % PROGRESS_EVERY == 0:
                log.info(
                    "%d reads processed (%d with MM), %.0f s elapsed",
                    stats.reads, stats.reads_with_mm, time.monotonic() - start,
                )
        if bar is not None:
            bar.close()
    if index_output:
        # Output is closed at this point, so the index matches the finished file.
        log.info("indexing %s", out_path)
        pysam.index(out_path)
    log.info(
        "done: %d reads (%d with MM) in %.1f s",
        stats.reads, stats.reads_with_mm, time.monotonic() - start,
    )
    return stats
