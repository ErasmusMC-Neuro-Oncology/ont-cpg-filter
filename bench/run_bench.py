#!/usr/bin/env python3
"""Time the filter with and without the CpG index, on the same BAM.

Variants:
    baseline    reads the BAM and fetches the tags, no filtering
    indexed     CpG bitmaps loaded from <fasta>.<MOTIF>_idx
    on-the-fly  CpG bitmaps built per chromosome during the run

Usage (from the project root):
    python bench/run_bench.py test_region.bam -r genome.fa
    python bench/run_bench.py test_region.bam -r genome.fa --variants indexed on-the-fly

The index is built before the timed runs if it is missing. If it exists but does
not match the FASTA or motif, the run stops with an error.
Timing covers reading, filtering and tag inspection, not writing a BAM.
Run it twice if the first run is slow: the second run has the OS file cache warm.
"""

from __future__ import annotations

import argparse
import contextlib
import functools
import os
import time

import pysam

from ont_cpg_filter.filter import Stats, filter_read
from ont_cpg_filter.index import build_index, index_path, load_index
from ont_cpg_filter.reference import Reference


def _no_filter(read, ref, stats):
    pass


VARIANTS = {
    "baseline": (lambda fasta, motif: contextlib.nullcontext(None), _no_filter),
    "indexed": (functools.partial(Reference, use_index=True), filter_read),
    "on-the-fly": (functools.partial(Reference, use_index=False), filter_read),
}


def run(make_ref, filt, bam, fasta, motif, limit):
    stats = Stats()
    out = []
    with make_ref(fasta, motif) as ref, pysam.AlignmentFile(bam, "rb") as inbam:
        used_index = getattr(ref, "uses_index", False)
        t0 = time.perf_counter()
        for i, read in enumerate(inbam.fetch(until_eof=True)):
            if limit and i >= limit:
                break
            stats.reads += 1
            if not read.is_unmapped and read.has_tag("MM"):
                stats.reads_with_mm += 1
                filt(read, ref, stats)
            mm = read.get_tag("MM") if read.has_tag("MM") else None
            ml = tuple(read.get_tag("ML")) if read.has_tag("ML") else None
            out.append((read.query_name, mm, ml))
        dt = time.perf_counter() - t0
    return dt, stats, out, used_index


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("bam")
    p.add_argument("-r", "--reference", required=True)
    p.add_argument("-m", "--motif", default="CG")
    p.add_argument("--limit", type=int, default=0, help="only the first N reads (0 = all)")
    p.add_argument("--variants", nargs="+", default=list(VARIANTS), choices=list(VARIANTS))
    a = p.parse_args()
    motif = a.motif.upper()

    # Preparation happens here, outside the timed runs.
    if not os.path.exists(a.reference + ".fai"):
        print("prep: creating FASTA index .fai")
        pysam.faidx(a.reference)
    if "indexed" in a.variants:
        existing = load_index(a.reference, motif)  # raises if stale or unreadable
        if existing is None:
            t0 = time.perf_counter()
            build_index(a.reference, motif)
            print(f"prep: built {index_path(a.reference, motif)} in {time.perf_counter() - t0:.1f} s (not timed)")
        else:
            existing.close()
            print(f"prep: existing {index_path(a.reference, motif)} matches FASTA and motif")

    results = {}
    for name in a.variants:
        make_ref, filt = VARIANTS[name]
        dt, st, out, used_index = run(make_ref, filt, a.bam, a.reference, motif, a.limit)
        results[name] = (dt, st, out)
        rate = st.reads / dt if dt else float("inf")
        note = ""
        if name == "indexed":
            note = " index=yes" if used_index else " index=NO"
        print(
            f"{name:11s} reads={st.reads} with_mm={st.reads_with_mm} "
            f"time={dt:.2f}s rate={rate:.0f} reads/s kept={st.kept} removed={st.removed} "
            f"added={st.added}{note}"
        )

    # Correctness: the indexed and on-the-fly variants must give identical output.
    if "indexed" in results and "on-the-fly" in results:
        a_out = results["indexed"][2]
        b_out = results["on-the-fly"][2]
        diff = [i for i, (x, y) in enumerate(zip(a_out, b_out)) if x != y]
        msg = f"on-the-fly vs indexed: {len(diff)} reads differ"
        if diff:
            msg += f", first: {a_out[diff[0]][0]}"
        print(msg)


if __name__ == "__main__":
    main()
