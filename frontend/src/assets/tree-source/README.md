# Tree source data

`generateTree.R` produces `tree.json` (nested) from NCBI taxonomy. The SPA
imports the flat `../taxa.json` and builds the hierarchy at runtime in
`src/components/tree/buildTree.ts`, so `tree.json` is not imported by any
code. It is kept here as the provenance of `taxa.json` and as the input to
any future server-side taxonomy endpoint. Do not delete it.
