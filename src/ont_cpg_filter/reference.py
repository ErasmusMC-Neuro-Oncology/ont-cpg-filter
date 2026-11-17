"""Reference access via pysam.FastaFile, with CpG bitmaps per chromosome.

The reference is the only source for CpG sites. For each chromosome two
bitmaps, packed to 1 bit per base, are used:
  - fwd: position is the first base of the motif (the C of CpG)
  - rev: position is the second base of the motif (the G of CpG)

The bitmaps come from <fasta>.<MOTIF>_idx if that file matches the FASTA and motif
(see index.py), otherwise they are built per chromosome when first needed.
"""

from __future__ import annotations

import logging
import os

import numpy as np
import pysam

from ont_cpg_filter.index import build_bitmaps, index_path, load_index

log = logging.getLogger(__name__)


class Reference:
    def __init__(self, fasta_path: str, motif: str = "CG", use_index: bool = True) -> None:
        if len(motif) != 2:
            raise ValueError("motif must be a dinucleotide, e.g. 'CG'")
        if not os.path.exists(fasta_path + ".fai"):
            log.info("no index found, creating %s.fai (this can take a while for a genome)", fasta_path)
            pysam.faidx(fasta_path)
        log.info("opening reference %s", fasta_path)
        self._fasta = pysam.FastaFile(fasta_path)
        self.motif = motif.upper()
        self._index = load_index(fasta_path, self.motif) if use_index else None
        if use_index and self._index is None:
            log.info(
                "no CpG index %s, computing bitmaps per chromosome",
                index_path(fasta_path, self.motif),
            )
        self._chrom: str | None = None
        self._len = 0
        self._fwd = np.zeros(0, dtype=np.uint8)
        self._rev = np.zeros(0, dtype=np.uint8)

    @property
    def uses_index(self) -> bool:
        """True if a matching <fasta>.<MOTIF>_idx was loaded."""
        return self._index is not None

    def close(self) -> None:
        if self._index is not None:
            self._index.close()
        self._fasta.close()

    def __enter__(self) -> "Reference":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def _load_chrom(self, chrom: str) -> None:
        """Make the bitmaps for one chromosome current (from the index or built now)."""
        if chrom == self._chrom:
            return
        if self._index is not None and chrom in self._index:
            self._len, self._fwd, self._rev = self._index.get(chrom)
        else:
            log.debug("building CpG bitmaps for %s", chrom)
            n, self._fwd, self._rev = build_bitmaps(self._fasta.fetch(chrom), self.motif)
            self._len = n
        self._chrom = chrom

    def site_mask(self, chrom: str, rpos: np.ndarray, forward: bool) -> np.ndarray:
        """Boolean mask: which reference positions in `rpos` are motif sites.

        forward=True: the position must be the first motif base (C of CpG).
        forward=False: the position must be the second motif base (G of CpG).
        A position of -1 (no reference position, e.g. soft clip) is never a site.
        """
        self._load_chrom(chrom)
        if self._len == 0:
            return np.zeros(rpos.shape, dtype=bool)
        packed = self._fwd if forward else self._rev
        ok = (rpos >= 0) & (rpos < self._len)
        idx = np.where(ok, rpos, 0)
        shift = (7 - (idx & 7)).astype(np.uint8)
        bits = (packed[idx >> 3] >> shift) & 1
        return ok & bits.astype(bool)

    def _fetch(self, chrom: str, start: int, end: int) -> str:
        """Uppercase reference sequence [start, end); empty if out of bounds."""
        if start < 0 or end > self._fasta.get_reference_length(chrom):
            return ""
        return self._fasta.fetch(chrom, start, end).upper()

    def keeps(self, chrom: str, pos: int, base: str) -> bool:
        """Scalar check, mainly for tests. True if a call on `base` at `pos` is in the motif.

        `base` is the base in the stored SEQ that carries the modification.
        """
        x, y = self.motif
        if base == x and self._fetch(chrom, pos, pos + 2) == self.motif:
            return True
        if base == y and pos >= 1 and self._fetch(chrom, pos - 1, pos + 1) == self.motif:
            return True
        return False
