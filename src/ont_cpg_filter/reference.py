"""Reference sequence access via pysam.FastaFile (htslib, uses the .fai index)."""

from __future__ import annotations

import logging
import os

import pysam

log = logging.getLogger(__name__)


class Reference:
    def __init__(self, fasta_path: str, motif: str = "CG") -> None:
        if len(motif) != 2:
            raise ValueError("motif must be a dinucleotide, e.g. 'CG'")
        if not os.path.exists(fasta_path + ".fai"):
            log.info("no index found, creating %s.fai (this can take a while for a genome)", fasta_path)
            pysam.faidx(fasta_path)
        log.info("opening reference %s", fasta_path)
        self._fasta = pysam.FastaFile(fasta_path)
        self.motif = motif.upper()

    def close(self) -> None:
        self._fasta.close()

    def __enter__(self) -> "Reference":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def _fetch(self, chrom: str, start: int, end: int) -> str:
        """Uppercase reference sequence [start, end); empty if out of bounds."""
        if start < 0 or end > self._fasta.get_reference_length(chrom):
            return ""
        return self._fasta.fetch(chrom, start, end).upper()

    def window(self, chrom: str, start: int, end: int) -> str:
        """Uppercase reference sequence [start, end) in one fetch (no bounds check)."""
        return self._fasta.fetch(chrom, start, end).upper()

    def keeps(self, chrom: str, pos: int, base: str) -> bool:
        """True if a modification call on `base` at 0-based reference `pos` is in the motif.

        `base` is the base in the stored SEQ that carries the modification
        ('C' for a call on the + strand, 'G' for a call on the - strand).
        If `base` is the first motif letter, the motif must start at `pos`
        (forward CpG). If it is the second letter, the motif must end at `pos`
        (reverse-strand call: the modified G is the G of the CpG).
        """
        x, y = self.motif
        if base == x and self._fetch(chrom, pos, pos + 2) == self.motif:
            return True
        if base == y and pos >= 1 and self._fetch(chrom, pos - 1, pos + 1) == self.motif:
            return True
        return False
