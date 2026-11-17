# ont-cpg-filter

Removes 5mC/5hmC modification calls from Oxford Nanopore BAM files at
positions that are not part of a CpG dinucleotide in the reference.

Modification calls live in the optional `MM` (which bases, which positions)
and `ML` (probabilities) tags. This tool rewrites those tags; the SAM flags
are not touched.

## Status

Skeleton only. CLI, reference access and tests are in place. MM/ML rewriting
is not implemented yet (see `filter.py`).

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

## Usage

```bash
ont-cpg-filter input.bam output.bam --reference genome.fa --motif CG
```

The reference FASTA must match the genome build the BAM was aligned to. An
index (`genome.fa.fai`) is created automatically if it is missing.

## Tests

```bash
pytest
```
