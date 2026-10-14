# intervis — project analysis records

Supporting analysis, figures and records for the SCIE3220 project
*"Development of a bioinformatics pipeline for identifying and visualising
integron cassette arrays in Vibrio genomes"* (intervis). The intervis tool
itself (engine, CLI, R/Shiny app) is in the repository root; this `records/`
tree holds the detection benchmark, the synthetic benchmark and the project
record keeping.

## Contents
- `record_keeping.md` — chronological research record for the project.
- `analysis/scripts/` — detection/scoring and synthetic-genome code:
  - `score_ap.py` — average-precision / precision–recall scoring (real genomes).
  - `make_synthetic.py` — synthetic-genome generator (known attC truth).
  - `run_synth_tools.sh` — runs IntegronFinder, HattCI and BacAnt on the synthetic genomes.
  - `score_synth.py` — scores detector output against the known synthetic coordinates.
  - (F1 scoring, `score_f1.py`, uses the same 100-bp midpoint matching and runs in the
    benchmark directory on Bunya.)
- `analysis/figures/` — evaluation figures and the scripts that make them
  (self-contained; matplotlib + numpy):
  - `fig2_detection_counts.*` — attC detected per tool, with reported-cassette reference lines.
  - `fig3_synthetic_recall.*` — recall vs attC divergence, random vs compensatory substitution.
  - `figA1_pr_curves.*` — HattCI precision–recall curves (average precision).
- `analysis/data/SCIE3220_detection_benchmark.xlsx` — benchmark workbook
  (reference genomes, detection counts, F1, average precision, synthetic benchmark).
- `report/` — current research report.

## Reproducing
See the "Reproducibility" section of `record_keeping.md` for tool versions,
exact commands, database release and accessions. Detection ran on UQ Bunya in
micromamba environments; the figure scripts run anywhere with matplotlib + numpy.
