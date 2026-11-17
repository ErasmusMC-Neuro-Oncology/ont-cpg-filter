"""Command-line interface."""

from __future__ import annotations

import argparse
import logging
import sys

from ont_cpg_filter import __version__
from ont_cpg_filter.filter import filter_bam


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ont-cpg-filter",
        description=(
            "Remove 5mC/5hmC modification calls (MM/ML tags) at positions "
            "that are not part of a CpG dinucleotide in the reference."
        ),
    )
    parser.add_argument("input", help="Input BAM (aligned, with MM/ML tags)")
    parser.add_argument("output", help="Output BAM")
    parser.add_argument(
        "-r",
        "--reference",
        required=True,
        help="Reference FASTA (same build as the BAM); .fai index is created if missing",
    )
    parser.add_argument(
        "--motif",
        default="CG",
        help="Dinucleotide that must contain the modified C (default: CG)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Threads for BAM (de)compression (default: 1)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="More logging (-v: debug)",
    )
    parser.add_argument(
        "--progress",
        action="store_true",
        help="show a progress bar (costs speed; off by default)",
    )
    parser.add_argument(
        "--no-output-index",
        action="store_true",
        help="do not write a .bai index for the output BAM",
    )
    parser.add_argument("--version", action="version", version=__version__)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    level = logging.DEBUG if args.verbose > 1 else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )
    try:
        stats = filter_bam(
            args.input,
            args.output,
            args.reference,
            motif=args.motif.upper(),
            threads=args.threads,
            index_output=not args.no_output_index,
            show_progress=args.progress,
        )
    except (NotImplementedError, ValueError, OSError) as exc:
        logging.getLogger("ont_cpg_filter").error("%s", exc)
        return 2
    logging.getLogger("ont_cpg_filter").info(
        "calls kept: %d, removed: %d", stats.kept, stats.removed
    )
    return 0
