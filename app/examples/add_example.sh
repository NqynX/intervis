#!/usr/bin/env bash
# Register one of your own genomes as a bundled intervis example.
#
#   ./add_example.sh <id> <genome.fna> [annotation.gff|.gbff]
#
# Copies the FASTA (and annotation, if given) into this examples/ folder and runs
# IntegronFinder once to precompute the .integrons, so the example loads instantly.
# The app auto-discovers examples in this folder at launch — the menu label is taken
# from the genome's FASTA header. Run inside the integronfinder env; relaunch the app.
#
# e.g.  ./add_example.sh mo6 ../../MO6.fasta ../../MO6.gb
set -euo pipefail
if [ "$#" -lt 2 ]; then
  echo "usage: $0 <id> <genome.fna> [annotation.gff|.gbff]"; exit 1
fi
ID="$1"; FNA="$2"; ANN="${3:-}"
DIR="$(cd "$(dirname "$0")" && pwd)"

cp "$FNA" "$DIR/${ID}.fna"
if [ -n "$ANN" ]; then EXT="${ANN##*.}"; cp "$ANN" "$DIR/${ID}.${EXT}"; fi

echo "[add_example] running IntegronFinder on ${ID}.fna (this can take a few minutes)..."
TMP="$(mktemp -d)"
integron_finder "$DIR/${ID}.fna" --outdir "$TMP" --local-max --keep-palindromes --cpu "${CPU:-8}"
cp "$(find "$TMP" -name '*.integrons' | head -1)" "$DIR/${ID}.integrons"
rm -rf "$TMP"

echo "[add_example] added '$ID'. Stop and relaunch the app — it appears in the Examples menu."
