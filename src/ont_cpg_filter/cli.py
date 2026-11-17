"""Command-line interface: `ont-cpg-filter filter ...` and `ont-cpg-filter index ...`."""

from __future__ import annotations

import argparse
import logging
import os
import sys
import time

from ont_cpg_filter import __version__
from ont_cpg_filter.filter import filter_bam
from ont_cpg_filter.index import build_index

log = logging.getLogger("ont_cpg_filter")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ont-cpg-filter",
        description=(
            "Keep only the base modification calls on CpG sites of the reference "
            "in the MM/ML tags of an aligned Oxford Nanopore BAM."
        ),
    )
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True, metavar="COMMAND")

    f = sub.add_parser(
        "filter",
        help="rewrite the MM/ML tags so that only CpG calls remain",
    )
    f.add_argument("input", help="input BAM (aligned, with MM/ML tags)")
    f.add_argument("output", help="output BAM")
    f.add_argument(
        "-r", "--reference", required=True,
        help="reference FASTA, same genome build as the BAM",
    )
    f.add_argument(
        "-m", "--motif", default="CG",
        help="two-base motif that must contain the modified base (default: CG)",
    )
    f.add_argument(
        "-t", "--threads", type=int, default=1,
        help="htslib threads for BAM compression and decompression (default: 1)",
    )
    f.add_argument(
        "-p", "--progress", action="store_true",
        help="show a progress bar (costs some speed)",
    )
    f.add_argument(
        "--no-output-index", action="store_true",
        help="do not write an index for the output BAM",
    )
    f.add_argument(
        "--no-reference-index", action="store_true",
        help="do not use <reference>.<MOTIF>_idx; compute the CpG bitmaps during the run",
    )
    f.add_argument("-v", "--verbose", action="count", default=0, help="debug logging with -vv")

    i = sub.add_parser(
        "index",
        help="build the CpG index <fasta>.<MOTIF>_idx for a reference",
    )
    i.add_argument("fasta", help="reference FASTA")
    i.add_argument("-m", "--motif", default="CG", help="two-base motif (default: CG)")
    i.add_argument("-v", "--verbose", action="count", default=0, help="debug logging with -vv")
    return parser


def _setup_logging(verbose: int) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose > 1 else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )


def _run_filter(args: argparse.Namespace, motif: str) -> int:
    stats = filter_bam(
        args.input,
        args.output,
        args.reference,
        motif=motif,
        threads=args.threads,
        index_output=not args.no_output_index,
        show_progress=args.progress,
        use_index=not args.no_reference_index,
    )
    log.info("calls kept: %d, removed: %d, added: %d", stats.kept, stats.removed, stats.added)
    return 0


def _run_index(args: argparse.Namespace, motif: str) -> int:
    t0 = time.perf_counter()
    out = build_index(args.fasta, motif)
    log.info(
        "built %s in %.1f s, %.1f MB",
        out, time.perf_counter() - t0, os.path.getsize(out) / 1e6,
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    _setup_logging(args.verbose)
    motif = args.motif.upper()
    if len(motif) != 2:
        log.error("motif must be a dinucleotide, e.g. CG")
        return 2
    try:
        if args.command == "index":
            return _run_index(args, motif)
        return _run_filter(args, motif)
    except (NotImplementedError, ValueError, OSError) as exc:
        log.error("%s", exc)
        return 2
