# Tree source data

`generateTree.R` writes `../taxa.json` (one NCBI lineage per species on
api.refgenie.org). Run it from this directory. The SPA builds the tree from
`taxa.json` at runtime in `src/components/tree/buildTree.ts`.
