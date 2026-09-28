# intervis

Interactive integron array visualiser. intervis takes a genome (it runs IntegronFinder
for you) or an existing IntegronFinder `.integrons` file and produces an interactive HTML
figure of the cassette array: pan and zoom, cassette arrows coloured by function, attC
sites, gene-name labels, and, for a two-genome comparison, ribbons that link shared
cassettes and show their percent identity.

It was written for the SCIE3220 project on integron cassette arrays in *Vibrio* and other
genomes.

<p align="center">
  <img src="docs/compare_view.png" alt="intervis comparison view" width="88%"><br>
  <em>Comparison view: two arrays aligned at the integrase, with ribbons joining homologous cassettes.</em>
</p>
<p align="center">
  <img src="docs/single_view.png" alt="intervis single-array view" width="88%"><br>
  <em>Single-array view: forward and reverse strands, attC sites, function colours, gene labels.</em>
</p>

## Features

- Runs IntegronFinder on a genome, or reads an existing `.integrons` file.
- Single-array view: forward and reverse strands, attC sites, cassettes coloured by
  function, and gene-name labels.
- Comparison view: two genomes aligned at the integrase, with ribbons joining homologous
  cassettes and showing percent identity on hover.
- Cassette naming from Swiss-Prot with DIAMOND, plus optional gene symbols from a GFF or
  GenBank file.
- Output is a single HTML file that runs in any browser with no server, and can be
  exported to SVG, PNG or PDF.
- A command-line tool for scripting and a Shiny web app for interactive use, both driven
  by the same engine.

## How it works

1. `intervis run` calls IntegronFinder (`--local-max --keep-palindromes`) on the genome
   and reads the `.integrons` file it writes.
2. The array is parsed, oriented with the integrase on the left, and each cassette is
   classified (AMR, toxin-antitoxin, anti-phage defence, other, hypothetical).
3. Cassettes are named with DIAMOND against Swiss-Prot. If a matching GFF or GenBank file
   is given, gene symbols are added by coordinate.
4. For one genome, intervis writes the single-array view. For two genomes, it clusters
   cassettes across them with DIAMOND reciprocal best hits and writes the comparison view
   with ribbons.

`intervis view` and `intervis compare` read `.integrons` files directly and skip
IntegronFinder, so they are fast. IntegronFinder is the only slow step, and only
`intervis run` calls it.

## Requirements

- Python 3.9 or newer (the engine uses the standard library only)
- [IntegronFinder](https://github.com/gem-pasteur/Integron_Finder) 2.0.2 or newer on PATH, to go from a genome to `.integrons`
- [DIAMOND](https://github.com/bbuchfink/diamond) 2.1 or newer on PATH, for cassette naming and comparison ribbons
- R 4.1 or newer with the `shiny` package, for the web app

## Installation

```bash
# 1. create the environment with IntegronFinder, DIAMOND and Python
micromamba create -f environment.yml      # or: conda env create -f environment.yml
micromamba activate integronfinder

# 2. install the intervis command-line tool
pip install -e .
```

To name cassettes, build a Swiss-Prot DIAMOND database once and pass it with
`--swissprot` (or set `INTERVIS_SWISSPROT`):

```bash
diamond makedb --in uniprot_sprot.fasta -d db/swissprot
```

## Command-line usage

```bash
# run IntegronFinder on a genome, then render the figure
intervis run genome.fna -o genome.html --cpu 8 --annotate --swissprot db/swissprot

# render directly from an existing .integrons file
intervis view N16961.integrons --genome N16961.fna -o N16961.html

# compare two genomes (the web app allows two; the command line allows more)
intervis compare A.integrons B.integrons --genomes A.fna B.fna -o comparison.html
```

To add gene names from a GFF or GenBank file, use `--annotation genome.gff` for a single
genome or `--annotations A.gff B.gff` for a comparison. A matching `.gff` or `.gbff` next
to the FASTA is found automatically.

## Web app

```bash
micromamba activate integronfinder
R -e 'shiny::runApp("app", host="0.0.0.0", port=8787, launch.browser=FALSE)'
```

Open the forwarded port in a browser. Pick a bundled example to see a result straight
away, or upload a genome (FASTA), optionally with a GFF or GenBank file for gene names,
and press Run. A `.integrons` file renders at once; a genome runs the full IntegronFinder
pipeline, which takes a few minutes.

Settings are read from environment variables, or from the defaults at the top of
`app/app.R`: `INTERVIS_PKG`, `INTERVIS_SWISSPROT`, `INTERVIS_CPU`.

## Examples

`app/examples/` includes two small synthetic examples (`DemoVibrioA` and `DemoVibrioB`)
that load instantly and show both views. Their comparison has five ribbons, one below
100 percent identity so the identity gradient is visible. The app finds any example in
this folder at startup, so examples you add stay available after updates.

Two helper scripts (run inside the environment):

```bash
# download real Vibrionaceae superintegron genomes from NCBI and build them as examples
cd app/examples && ./fetch_examples.sh

# add one of your own genomes as an example (precomputes its .integrons)
./add_example.sh <id> <genome.fna> [annotation.gff|.gbff]
```

`fetch_examples.sh` keeps a genome only if IntegronFinder finds a complete integron (an
integrase and a real attC array), so every example is comparable. Downloaded genomes are
ignored by git; only the small examples are committed.

## Running on an HPC (Bunya)

The app needs IntegronFinder, DIAMOND and Python, so it runs on a compute node rather than
a hosted Shiny service. Request an interactive allocation, start the app on the node, and
forward the port:

```bash
# on the login node
salloc --account=<account> --partition=general --nodes=1 --ntasks=1 \
       --cpus-per-task=8 --mem=32G --time=04:00:00 --job-name=intervis
srun --export=ALL --pty bash -l
hostname          # note the node, for example bun137

micromamba activate integronfinder
R -e 'shiny::runApp("app", host="0.0.0.0", port=8787, launch.browser=FALSE)'
```

```bash
# from your laptop, using the node name from above
ssh -N -L 8787:bun137:8787 <user>@bunya.rcc.uq.edu.au
# then open http://localhost:8787
```

## Repository layout

```
intervis/
  intervis/            Python engine and CLI (standard library only)
    cli.py             run | view | compare
    engine.py          IntegronFinder wrapper, .integrons parser, HTML builder
    annotate.py        cassette naming (DIAMOND/Swiss-Prot; GFF/GenBank symbols)
    cluster.py         cross-genome cassette clustering (DIAMOND reciprocal best hits)
    templates/         cartographer.html (single), synteny.html (comparison)
  app/                 R/Shiny web app
    app.R
    examples/          synthetic examples, fetch_examples.sh, add_example.sh
  environment.yml      conda/micromamba environment
  pyproject.toml
  docs/                screenshots
```

## Status and plans

Working now: single and comparison views, cassette naming (Swiss-Prot and optional
GFF/GenBank symbols), cross-genome ribbons with percent identity, bundled examples, and
the Shiny app.

Planned: package the app as an installable R Shiny package, and add population-scale
cassette context (LexicMap / AllTheBacteria) as a command-line feature.

## Acknowledgements

intervis uses [IntegronFinder](https://github.com/gem-pasteur/Integron_Finder) (Néron et
al.) to find integrons and [DIAMOND](https://github.com/bbuchfink/diamond) (Buchfink et
al.) to cluster cassette proteins. The comparison view is modelled on
[clinker](https://github.com/gamcil/clinker) (Gilchrist and Chooi).

## License

MIT. See [LICENSE](LICENSE).

## Citation

If you use intervis, please cite it along with IntegronFinder and DIAMOND:

> Nguyen, L. (2026). intervis: an interactive integron array visualiser. SCIE3220.
> https://github.com/NqynX/intervis
