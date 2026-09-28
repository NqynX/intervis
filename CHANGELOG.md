# Changelog

All notable changes to intervis. This project loosely follows
[Semantic Versioning](https://semver.org/).

## [1.0.0] — 2026-09-28

First public release. A complete integron array visualiser: engine + CLI with a thin
R/Shiny GUI.

### Engine / CLI
- `intervis run` — genome → IntegronFinder (`--local-max --keep-palindromes`) → parse →
  annotate → interactive HTML (single genome, or a comparison for two).
- `intervis view` / `intervis compare` — render instantly from existing `.integrons`
  files, no IntegronFinder needed.
- Cassette naming via DIAMOND vs Swiss-Prot; optional real gene symbols from a matching
  GFF3/GenBank overlaid by coordinate. Cassettes classified as AMR / toxin–antitoxin /
  anti-phage defence / other / hypothetical.
- Cross-genome cassette clustering by DIAMOND **reciprocal best hits**, so a genome vs
  itself gives 100 % identity for every cassette's twin and ribbons carry exact per-pair
  identity + coverage.
- Base-level nucleotide zoom (single and comparison views).

### Single-array viewer (`cartographer.html`)
- Two-strand layout (forward above / reverse below the backbone), attC sites as
  toggleable asterisks, function colours, clinker-style gene-name **Labels** (italic only
  for true gene symbols), per-arrow recolour, region shading, pan/zoom, and a left-gutter
  genome identifier.
- Export to SVG / PNG / PDF (no baked-in title; genome name in the filename).

### Comparison viewer (`synteny.html`)
- Same look as the single view, stacked, with two genomes aligned at the integrase and
  homologous cassettes joined by filled ribbons (coloured by homology group or by exact
  identity, toggleable). Ribbons attach to the arrows, are clipped under zoom, and show
  exact % identity on hover.

### GUI (`app/app.R`)
- Upload 1–2 genomes (FASTA) with optional GFF/GenBank, or a `.integrons` for an instant
  render; upload cap raised to 500 MB.
- Bundled examples auto-discovered from `app/examples/` (survive upgrades); ships two
  synthetic demos, with `fetch_examples.sh` (real Vibrionaceae superintegrons) and
  `add_example.sh` (register your own genome) as helpers. `fetch_examples.sh` keeps a
  genome only if IntegronFinder finds a complete integron (integrase + attC array).

[1.0.0]: https://github.com/NqynX/intervis/releases/tag/v1.0.0
