#!/bin/bash

set -euo pipefail

PY=$(compgen -c python3. | grep -E '^python3\.[0-9]+$' | sort -uV | tail -n1)
PY=${PY:-python3}   # val terug op python3 als er geen python3.X gevonden is
echo "Gebruik $PY ($($PY --version))"

"$PY" -m venv .venv

source .venv/bin/activate
pip install -e ".[dev]"
pytest
