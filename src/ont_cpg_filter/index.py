"""Serialized CpG bitmap index: <fasta>.cg_idx

Header checks are deliberately minimal: FASTA size and mtime, plus the motif
(the bitmaps depend on it). No checksums and no contig name check.

Build once:
    python -m ont_cpg_filter.index genome.fa --motif CG

Memory: build_index keeps all bitmaps in memory before writing. For a human
genome that is roughly 2 x 3.1 Gb / 8 bytes, about 780 MB.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import time

import numpy as np
import pysam

log = logging.getLogger(__name__)

VERSION = 1


def index_path(fasta_path: str) -> str:
    return fasta_path + ".cg_idx"


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
    out = out_path or index_path(fasta_path)
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
    """Lazy access to the bitmaps in a .cg_idx file (one contig loaded at a time)."""

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
    """Return the index, or None if there is no index file.

    Raises ValueError if the file exists but is unreadable or does not match
    the FASTA size, mtime and motif. A stale index is a problem to fix, not a
    reason to silently fall back.
    """
    path = index_path(fasta_path)
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
        or header.get("motif") != motif
        or header.get("fasta_size") != stamp["fasta_size"]
        or header.get("fasta_mtime_ns") != stamp["fasta_mtime_ns"]
    ):
        z.close()
        raise ValueError(
            f"{path} is stale or was built for another motif (header motif "
            f"{header.get('motif')!r}, requested {motif!r}). Rebuild it with: "
            f"python -m ont_cpg_filter.index {fasta_path} --motif {motif}"
        )
    log.info("using CpG index %s", path)
    return CpgIndex(z, header["contigs"])


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Build the <fasta>.cg_idx CpG index.")
    p.add_argument("fasta")
    p.add_argument("--motif", default="CG")
    p.add_argument("-o", "--output", help="default: <fasta>.cg_idx")
    a = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    t0 = time.perf_counter()
    out = build_index(a.fasta, a.motif.upper(), a.output)
    log.info("built %s in %.1f s, %.1f MB", out, time.perf_counter() - t0, os.path.getsize(out) / 1e6)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
