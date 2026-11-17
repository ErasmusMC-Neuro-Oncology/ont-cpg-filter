"""Version 2: one reference fetch per read, CpG sites as Python sets from a regex scan."""

from __future__ import annotations

import re
from array import array

import pysam

from ont_cpg_filter.filter import parse_mm, serialize_mm

_COMP = str.maketrans("ACGTN", "TGCAN")


class FastaWindow:
    def __init__(self, fasta_path: str, motif: str = "CG") -> None:
        self.f = pysam.FastaFile(fasta_path)
        self.motif = motif.upper()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.f.close()

    def window(self, chrom: str, start: int, end: int) -> str:
        return self.f.fetch(chrom, start, end).upper()


def _sites(refseq: str, motif: str, start: int) -> tuple[set[int], set[int]]:
    starts = [m.start() for m in re.finditer(f"(?={motif})", refseq)]
    return {start + p for p in starts}, {start + p + 1 for p in starts}


def filter_read(read: pysam.AlignedSegment, ref: FastaWindow, stats) -> None:
    seq = read.query_sequence
    if seq is None:
        return
    entries = parse_mm(read.get_tag("MM"))
    ml_in = list(read.get_tag("ML")) if read.has_tag("ML") else None
    chrom = read.reference_name
    # Extend the window by one base on each side: a CpG at the read edge needs
    # the neighbouring base to be in the window (C at rs-1 for a G at rs, G at re for a C at re-1).
    ws = max(read.reference_start - 1, 0)
    we = min(read.reference_end + 1, ref.f.get_reference_length(chrom))
    fwd, rev = _sites(ref.window(chrom, ws, we), ref.motif, ws)
    qpos_to_rpos = read.get_reference_positions(full_length=True)

    new_entries = []
    ml_out = []
    ml_pos = 0
    for entry in entries:
        strand = entry.strand
        if read.is_reverse:
            strand = "-" if strand == "+" else "+"
        seq_base = entry.base if strand == "+" else entry.base.translate(_COMP)
        cand = [m.start() for m in re.finditer(seq_base, seq)]
        if strand == "-":
            cand.reverse()
        n_codes = len(entry.codes)

        calls = []
        k = 0
        for s in entry.skips:
            k += s
            calls.append(k)
            k += 1
        if calls and calls[-1] >= len(cand):
            raise ValueError(f"{read.query_name}: MM skip counts exceed candidates")

        n_vals = len(calls) * n_codes
        vals = ml_in[ml_pos : ml_pos + n_vals] if ml_in is not None else None
        ml_pos += n_vals

        sites = fwd if seq_base == ref.motif[0] else rev
        kept = []
        for j, ci in enumerate(calls):
            rpos = qpos_to_rpos[cand[ci]]
            if rpos is not None and rpos in sites:
                kept.append(j)

        removed = len(calls) - len(kept)
        stats.kept += len(kept) * n_codes
        stats.removed += removed * n_codes
        if not kept:
            continue

        new_skips = []
        prev = -1
        for j in kept:
            new_skips.append(calls[j] - prev - 1)
            prev = calls[j]
        entry.skips = new_skips
        if removed:
            entry.flag = "?"
        new_entries.append(entry)
        if vals is not None:
            for j in kept:
                ml_out.extend(vals[j * n_codes : (j + 1) * n_codes])

    if new_entries:
        read.set_tag("MM", serialize_mm(new_entries))
        if ml_in is not None:
            read.set_tag("ML", array("B", ml_out))
    else:
        read.set_tag("MM", None)
        if read.has_tag("ML"):
            read.set_tag("ML", None)
