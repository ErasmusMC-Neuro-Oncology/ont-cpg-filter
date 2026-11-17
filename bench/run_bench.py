#!/usr/bin/env python3
"""Run the filter implementations on the same BAM, time them and compare outputs.

Versions:
    v1-percall   one pysam fetch per CpG check
    v2-window    one reference fetch per read, regex scan for sites
    v3-index     numpy bitmaps loaded from <fasta>.cg_idx
    v3-noindex   numpy bitmaps built per chromosoom during the run

Usage (from the project root):
    python bench/run_bench.py test_region.bam -r genome.fa
    python bench/run_bench.py test_region.bam -r genome.fa --versions v3-index v3-noindex

Build the index first for v3-index: python -m ont_cpg_filter.index genome.fa
Timing covers reading, filtering and tag inspection, not writing a BAM.
Run it twice if the first run is slow: the second run has the OS file cache warm.
"""

from __future__ import annotations

import argparse
import functools
import os
import time

import pysam

import v1_percall
import v2_window
import v3_bitmap
from ont_cpg_filter.filter import Stats
from ont_cpg_filter.index import build_index, load_index

VERSIONS = {
    "v1-percall": (v1_percall.FastaKeeps, v1_percall.filter_read),
    "v2-window": (v2_window.FastaWindow, v2_window.filter_read),
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
        existing = load_index(a.reference, motif)
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
        dt, st, out, used_index = run(make_ref, filt, a.bam, a.reference, a.motif.upper(), a.limit)
        results[name] = (dt, st, out)
        rate = st.reads / dt if dt else float("inf")
        note = ""
        if name == "v3-index":
            note = " index=yes" if used_index else " index=NO (no matching .cg_idx, bitmaps were built)"
        print(
            f"{name:12s} reads={st.reads} with_mm={st.reads_with_mm} "
            f"time={dt:.2f}s rate={rate:.0f} reads/s kept={st.kept} removed={st.removed}{note}"
        )

    # v1 is the most direct implementation of the rule, so it is the reference for correctness.
    base = "v1-percall" if "v1-percall" in results else a.versions[0]
    base_out = results[base][2]
    for name in a.versions:
        if name == base:
            continue
        diff = [i for i, (x, y) in enumerate(zip(base_out, results[name][2])) if x != y]
        msg = f"{name} vs {base}: {len(diff)} reads differ"
        if diff:
            msg += f", first: {base_out[diff[0]][0]}"
        print(msg)


if __name__ == "__main__":
    main()
