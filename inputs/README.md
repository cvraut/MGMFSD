# Inputs

This directory contains a bunch of possible input files.

## Format 1 - xlsx
- 2 columns
  - 1st col is gene name
  - 2nd col is guide sequence
    - guide sequence is always assumed to be 5' --> 3'
      - primer (NGG or G) will always try to match on the 3' side (so at the end of the guide strand)
- N rows
  - N (-1) genes
  - first row may be header
- for multiple guides in same gene, put unique guides on unique rows