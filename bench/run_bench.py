#!/usr/bin/env python3
"""Time the v3 filter with and without the CpG index, on the same BAM.

Versions:
    v0-passthrough  reads the BAM and fetches the tags, no filtering (baseline)
    v3-index        bitmaps loaded from <fasta>.cg_idx
    v3-noindex      bitmaps built per chromosome during the run

Usage (from the project root):
    python bench/run_bench.py test_region.bam -r genome.fa
    python bench/run_bench.py test_region.bam -r genome.fa --versions v3-index v3-noindex

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

import v3_bitmap
from ont_cpg_filter.filter import Stats
from ont_cpg_filter.index import build_index, load_index


def _no_filter(read, ref, stats):
    pass


VERSIONS = {
    "v0-passthrough": (lambda fasta, motif: contextlib.nullcontext(None), _no_filter),
    "v3-index": (functools.partial(v3_bitmap.FastaRef, use_index=True), v3_bitmap.filter_read),
    "v3-noindex": (functools.partial(v3_bitmap.FastaRef, use_index=False), v3_bitmap.filter_read),
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
    p.add_argument("--motif", default="CG")
    p.add_argument("--limit", type=int, default=0, help="only the first N reads (0 = all)")
    p.add_argument("--versions", nargs="+", default=list(VERSIONS), choices=list(VERSIONS))
    a = p.parse_args()
    motif = a.motif.upper()

    # Preparation happens here, outside the timed runs.
    if not os.path.exists(a.reference + ".fai"):
        print("prep: creating FASTA index .fai")
        pysam.faidx(a.reference)
    if "v3-index" in a.versions:
        existing = load_index(a.reference, motif)  # raises if stale or unreadable
        if existing is None:
            t0 = time.perf_counter()
            build_index(a.reference, motif)
            print(f"prep: built {a.reference}.cg_idx in {time.perf_counter() - t0:.1f} s (not timed)")
        else:
            existing.close()
            print(f"prep: existing {a.reference}.cg_idx matches FASTA and motif")

    results = {}
    for name in a.versions:
        make_ref, filt = VERSIONS[name]
        dt, st, out, used_index = run(make_ref, filt, a.bam, a.reference, motif, a.limit)
        results[name] = (dt, st, out)
        rate = st.reads / dt if dt else float("inf")
        note = ""
        if name == "v3-index":
            note = " index=yes" if used_index else " index=NO"
        print(
            f"{name:15s} reads={st.reads} with_mm={st.reads_with_mm} "
            f"time={dt:.2f}s rate={rate:.0f} reads/s kept={st.kept} removed={st.removed} "
            f"added={st.added}{note}"
        )

    # Correctness: the two v3 variants must give identical output.
    if "v3-index" in results and "v3-noindex" in results:
        a_out = results["v3-index"][2]
        b_out = results["v3-noindex"][2]
        diff = [i for i, (x, y) in enumerate(zip(a_out, b_out)) if x != y]
        msg = f"v3-noindex vs v3-index: {len(diff)} reads differ"
        if diff:
            msg += f", first: {a_out[diff[0]][0]}"
        print(msg)


if __name__ == "__main__":
    main()
