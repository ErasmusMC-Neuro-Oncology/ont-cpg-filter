"""Version 1: one pysam fetch per CpG check (strand fix included)."""

from __future__ import annotations

from array import array

import pysam

from ont_cpg_filter.filter import parse_mm, serialize_mm

_COMP = str.maketrans("ACGTN", "TGCAN")


class FastaKeeps:
    def __init__(self, fasta_path: str, motif: str = "CG") -> None:
        self.f = pysam.FastaFile(fasta_path)
        self.motif = motif.upper()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.f.close()

    def _fetch(self, chrom: str, start: int, end: int) -> str:
        if start < 0 or end > self.f.get_reference_length(chrom):
            return ""
        return self.f.fetch(chrom, start, end).upper()

    def keeps(self, chrom: str, pos: int, base: str) -> bool:
        x, y = self.motif
        if base == x and self._fetch(chrom, pos, pos + 2) == self.motif:
            return True
        if base == y and pos >= 1 and self._fetch(chrom, pos - 1, pos + 1) == self.motif:
            return True
        return False


def filter_read(read: pysam.AlignedSegment, ref: FastaKeeps, stats) -> None:
    seq = read.query_sequence
    if seq is None:
        return
    entries = parse_mm(read.get_tag("MM"))
    ml_in = list(read.get_tag("ML")) if read.has_tag("ML") else None
    qpos_to_rpos = {q: r for q, r in read.get_aligned_pairs(matches_only=True)}
    chrom = read.reference_name

    new_entries = []
    ml_out = []
    ml_pos = 0
    for entry in entries:
        strand = entry.strand
        if read.is_reverse:
            strand = "-" if strand == "+" else "+"
        seq_base = entry.base if strand == "+" else entry.base.translate(_COMP)
        cand = [i for i, b in enumerate(seq) if b == seq_base]
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

        kept = []
        for j, ci in enumerate(calls):
            rpos = qpos_to_rpos.get(cand[ci])
            if rpos is not None and ref.keeps(chrom, rpos, seq_base):
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
