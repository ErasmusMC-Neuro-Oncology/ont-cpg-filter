"""Serialized CpG bitmap index: <fasta>.<MOTIF>_idx, for example genome.fa.CG_idx

Header checks are deliberately minimal: FASTA size and mtime, plus the motif
(the bitmaps depend on it). No checksums and no contig name check.

Memory: build_index keeps all bitmaps in memory before writing. For a human
genome that is roughly 2 x 3.1 Gb / 8 bytes, about 780 MB.
"""

from __future__ import annotations

import json
import logging
import os

import numpy as np
import pysam

log = logging.getLogger(__name__)

VERSION = 1


def index_path(fasta_path: str, motif: str) -> str:
    return f"{fasta_path}.{motif.upper()}_idx"


def _fasta_stamp(fasta_path: str) -> dict:
    st = os.stat(fasta_path)  # follows symlinks, so the target's size and mtime are used
    return {"fasta_size": st.st_size, "fasta_mtime_ns": st.st_mtime_ns}


def build_bitmaps(seq: str, motif: str) -> tuple[int, np.ndarray, np.ndarray]:
    """Return (length, packed fwd bits, packed rev bits) for one contig."""
    arr = np.frombuffer(seq.upper().encode("ascii"), dtype=np.uint8)
    x, y = (ord(c) for c in motif)
    n = arr.size
    fwd = np.zeros(n, dtype=bool)
    rev = np.zeros(n, dtype=bool)
    hit = np.flatnonzero((arr[:-1] == x) & (arr[1:] == y))
    fwd[hit] = True
    rev[hit + 1] = True
    return n, np.packbits(fwd), np.packbits(rev)


def build_index(fasta_path: str, motif: str, out_path: str | None = None) -> str:
    motif = motif.upper()
    if len(motif) != 2:
        raise ValueError("motif must be a dinucleotide, e.g. 'CG'")
    out = out_path or index_path(fasta_path, motif)
    contigs = []
    arrays: dict[str, np.ndarray] = {}
    with pysam.FastaFile(fasta_path) as f:
        for i, chrom in enumerate(f.references):
            n, pf, pr = build_bitmaps(f.fetch(chrom), motif)
            contigs.append({"name": chrom, "length": n})
            arrays[f"fwd_{i}"] = pf
            arrays[f"rev_{i}"] = pr
            log.debug("built %s (%d bases)", chrom, n)
    header = {"version": VERSION, "motif": motif, "contigs": contigs, **_fasta_stamp(fasta_path)}
    arrays["header"] = np.frombuffer(json.dumps(header).encode(), dtype=np.uint8)
    tmp = out + ".tmp"
    with open(tmp, "wb") as fh:
        np.savez(fh, **arrays)
    os.replace(tmp, out)
    return out


class CpgIndex:
    """Lazy access to the bitmaps in an index file (one contig loaded at a time)."""

    def __init__(self, npz, contigs: list[dict]) -> None:
        self._z = npz
        self._pos = {c["name"]: (i, c["length"]) for i, c in enumerate(contigs)}

    def __contains__(self, chrom: str) -> bool:
        return chrom in self._pos

    def get(self, chrom: str) -> tuple[int, np.ndarray, np.ndarray]:
        i, n = self._pos[chrom]
        return n, self._z[f"fwd_{i}"], self._z[f"rev_{i}"]

    def close(self) -> None:
        self._z.close()


def load_index(fasta_path: str, motif: str) -> CpgIndex | None:
    """Return the index, or None if there is no index file for this motif.

    Raises ValueError if the file exists but is unreadable or does not match
    the FASTA size, mtime and motif. A stale index is a problem to fix, not a
    reason to silently fall back.
    """
    path = index_path(fasta_path, motif)
    if not os.path.exists(path):
        return None
    try:
        z = np.load(path)
        header = json.loads(bytes(z["header"]).decode())
    except Exception as exc:
        raise ValueError(f"cannot read {path}: {exc}") from exc
    stamp = _fasta_stamp(fasta_path)
    if (
        header.get("version") != VERSION
        or header.get("motif") != motif.upper()
        or header.get("fasta_size") != stamp["fasta_size"]
        or header.get("fasta_mtime_ns") != stamp["fasta_mtime_ns"]
    ):
        z.close()
        raise ValueError(
            f"{path} is stale or was built for another motif. "
            f"Rebuild it with: ont-cpg-filter index {fasta_path} --motif {motif.upper()}"
        )
    log.info("using CpG index %s", path)
    return CpgIndex(z, header["contigs"])
