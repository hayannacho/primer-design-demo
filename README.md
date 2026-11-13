# Primer Design Demo

Design a primer pair around a genomic position with a 600–800 bp product using Primer3.

## Results and verification

3 tests: primers bind the expected forward/reverse-complement sequence, product length and target coverage are correct, FASTA lookup works, and short windows/invalid positions are handled. The synthetic sequence yields one primer pair.

The bundled `example-output.json` records a local run. Examples are synthetic. Verified with Python 3.11.15 on macOS (CPU); floating-point results may vary slightly by platform.

## Run

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python model.py
python -m unittest discover -s tests -v
```

For a local FASTA file:

```bash
python model.py --fasta genome.fa --chromosome chr1 --position 1000
```

Coordinates in the output are one-based and inclusive. `pairs_found: 0` includes an explanation.

## Method and scope

The demo preserves the original one-based position, ±500 bp window and 200-base target convention. This edition uses primer3-py instead of subprocess/pyfasta and adds explicit no-pair reporting, validation and binding-coordinate checks. Defaults are a synthetic 2,000-base sequence. No specificity search or experimental validation is performed.

## Project context

BIOL 6505 · Fall 2025 · team project. The coursework integrates FASTA lookup, target-window selection, Primer3 configuration and primer extraction.
