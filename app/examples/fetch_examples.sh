#!/usr/bin/env bash
# Download real superintegron genomes from NCBI and build them as intervis examples
# (FASTA + GenBank + precomputed .integrons), for the single and compare views.
#
# Every example is KEPT ONLY IF IntegronFinder finds a COMPLETE integron — an integron
# INTEGRASE (intI) together with a real attC cassette array (>= MIN_ATTC sites). So the
# menu can never contain an integrase-less array (a CALIN, e.g. V. vulnificus CMCP6) or
# a tiny 1-2 cassette integron: all examples are comparable superintegrons with intI.
#
# Scope note: superintegrons of this scale (tens-to-hundreds of cassettes, ~50-140 kb)
# are a Vibrionaceae hallmark. This set spans that family across GENERA and SPECIES —
# Vibrio (cholerae, parahaemolyticus, vulnificus) and Photobacterium — rather than
# distant genera (Xanthomonas, Nitrosomonas, Shewanella), whose small sedentary /
# integrase-only chromosomal integrons are neither size-comparable nor always
# array-bearing, so they cannot be compared meaningfully.
#
# Run ON BUNYA, inside the integronfinder env (integron_finder + curl on PATH):
#     micromamba activate integronfinder
#     cd /scratch/user/uqcngu19/vibrio-integron-pipeline/intervis/app/examples
#     ./fetch_examples.sh
# It runs IntegronFinder on each genome, so budget a few minutes per genome.
set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
CPU="${CPU:-8}"
MIN_ATTC="${MIN_ATTC:-15}"   # keep only real arrays (superintegron-scale), not 1-2 cassette integrons

# id | NCBI accession(s). Comma-separate accessions to build ONE genome from several
# replicons (used where the superintegron's replicon is uncertain — IntegronFinder finds
# it on whichever). Menu labels are read from each genome's own FASTA header.
EXAMPLES=(
  "vcN16961|NC_002506.1"          # Vibrio cholerae N16961, chr II — the REFERENCE-STANDARD superintegron
  "vpRIMD|NC_004603.1"            # Vibrio parahaemolyticus RIMD 2210633, chr 1
  "vvYJ016|BA000037"             # Vibrio vulnificus YJ016, chr I
  "ppSS9|CR354531,CR354532"       # Photobacterium profundum SS9 (both chromosomes) — different genus, same family
)

command -v integron_finder >/dev/null || { echo "integron_finder not on PATH — run: micromamba activate integronfinder"; exit 1; }
command -v curl            >/dev/null || { echo "curl not found"; exit 1; }

# drop the non-comparable non-Vibrio attempts and the integrase-less V. vulnificus
# CMCP6 (IntegronFinder called it a CALIN — no integrase). The bundled synthetic
# demos (DemoVibrioA/B) are deliberately KEPT — they ship in the app for an instant,
# no-IntegronFinder demo of the single and compare views.
rm -f "$DIR"/xccampestris* "$DIR"/neuropaea* "$DIR"/shewMR1* \
      "$DIR"/vvCMCP6* "$DIR"/examples.tsv 2>/dev/null || true

# a genome is a usable example only if its .integrons has a COMPLETE integron
# (integrase + array) and enough attC sites to be superintegron-scale
complete_ok(){  # .integrons path -> 0 if a complete, large-enough integron is present
  local f="$1" attc
  [ -s "$f" ] || return 1
  grep -qw complete "$f" || return 1               # "complete" = intI + attC (vs CALIN=attC only, In0=intI only)
  attc="$(grep -c attC "$f" 2>/dev/null || echo 0)"
  [ "$attc" -ge "$MIN_ATTC" ]
}

fetch(){  # "acc1,acc2,..."  outbase   -> outbase.fna, outbase.gb (concatenated over replicons)
  local accs="$1" out="$2" base="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi" a
  : > "${out}.fna"; : > "${out}.gb"
  IFS=',' read -ra ACCS <<< "$accs"
  for a in "${ACCS[@]}"; do
    echo "[fetch] $a (FASTA + GenBank) ..."
    curl -fsS "${base}?db=nuccore&id=${a}&rettype=fasta&retmode=text"       >> "${out}.fna"
    curl -fsS "${base}?db=nuccore&id=${a}&rettype=gbwithparts&retmode=text" >> "${out}.gb"
  done
  [ -s "${out}.fna" ] || { echo "  download failed for $accs"; exit 1; }
}

for row in "${EXAMPLES[@]}"; do
  IFS='|' read -r id accs <<< "$row"
  if complete_ok "$DIR/$id.integrons"; then
    echo "[skip] $id already built (rm $id.* to rebuild)"; continue; fi
  fetch "$accs" "$DIR/$id"
  echo "[IntegronFinder] $id — this takes a few minutes ..."
  TMP="$(mktemp -d)"
  integron_finder "$DIR/$id.fna" --outdir "$TMP" --local-max --keep-palindromes --cpu "$CPU"
  cp "$(find "$TMP" -name '*.integrons' | head -1)" "$DIR/$id.integrons" 2>/dev/null || true
  rm -rf "$TMP"
  # keep ONLY a complete, comparable superintegron: integrase (intI) + a real attC array
  if complete_ok "$DIR/$id.integrons"; then
    echo "[done] $id  (complete superintegron, $(grep -c attC "$DIR/$id.integrons") attC sites)"
  else
    attc="$(grep -c attC "$DIR/$id.integrons" 2>/dev/null || echo 0)"
    echo "[drop] $id — need a COMPLETE integron (integrase intI + >= $MIN_ATTC attC sites); found attC=$attc, complete=$(grep -qw complete "$DIR/$id.integrons" 2>/dev/null && echo yes || echo no). Not a comparable superintegron — removing."
    rm -f "$DIR/$id".fna "$DIR/$id".gb "$DIR/$id".integrons
  fi
done

echo
echo "Built. Restart the app (Ctrl-C, relaunch) — every kept genome above appears in the"
echo "Examples menu. All carry an integron integrase + a real attC array, so any two are"
echo "directly comparable."
