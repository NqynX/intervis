#!/bin/bash
# =============================================================================
# run_synth_tools.sh  -  run the detectors on every synthetic genome
# Run from /scratch/user/uqcngu19/tool-benchmark/ AFTER make_synthetic.py.
# Submit as a batch job (several IntegronFinder runs take a while):
#     sbatch --account=a_uqccr --cpus-per-task=4 --mem=8G --time=06:00:00 \
#            --wrap "bash run_synth_tools.sh"
# =============================================================================
set -uo pipefail

# --- make micromamba usable inside a batch shell (it often isn't on PATH) ---
if ! command -v micromamba >/dev/null 2>&1; then
  for c in "$MAMBA_EXE" "$HOME/bin/micromamba" "$HOME/.local/bin/micromamba" \
           "$HOME/micromamba/bin/micromamba" "$HOME/micromamba"; do
    [ -n "${c:-}" ] && [ -x "$c" ] && export PATH="$(dirname "$c"):$PATH" && break
  done
fi
eval "$(micromamba shell hook --shell bash 2>/dev/null)" 2>/dev/null || true
if ! command -v micromamba >/dev/null 2>&1; then
  echo "ERROR: micromamba not found. Edit the path block at the top of this script." >&2
  exit 1
fi

SYN=synthetic          # directory holding synthetic_*.fna from make_synthetic.py
OUT=synth_runs
CPU=4
mkdir -p "$OUT"
shopt -s nullglob

n=0
for fna in "$SYN"/synthetic_*.fna; do
  g=$(basename "$fna" .fna); n=$((n+1))
  echo "==================== $g ===================="

  micromamba run -n integronfinder integron_finder \
      --circ --cpu "$CPU" --outdir "$OUT/IFdef_$g" "$fna" >/dev/null 2>&1 \
      && echo "  IF default  ok" || echo "  IF default  FAILED"

  micromamba run -n integronfinder integron_finder \
      --local-max --circ --cpu "$CPU" --outdir "$OUT/IFloc_$g" "$fna" >/dev/null 2>&1 \
      && echo "  IF local    ok" || echo "  IF local    FAILED"

  micromamba run -n hattci hattci.out -b "$fna" "$OUT/hattci_$g.txt" >/dev/null 2>&1 \
      && echo "  HattCI      ok" || echo "  HattCI      FAILED"

  micromamba run -n bacant bacant -n "$fna" -o "$OUT/bacant_$g" >/dev/null 2>&1 \
      && echo "  BacAnt      ok" || echo "  BacAnt      FAILED"
done
echo "processed $n genomes -> $OUT/"
