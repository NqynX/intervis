# SCIE3220 Research Project — Record Keeping

**Student:** Liam Nguyen (48596804)
**Project:** Development of a bioinformatics pipeline for identifying and visualising integron cassette arrays in *Vibrio* genomes (intervis)
**Supervisors:** Jan Engelstädter; Andrew Letten (Engelstädter Lab)
**Code and data:** development repository `vibrio-integron-pipeline`; released tool https://github.com/NqynX/intervis
**Compute:** UQ Bunya HPC cluster (micromamba environments, SLURM)

*A chronological record of a computational/bioinformatics project, grouped by work phase and dated to match the version-control history. Each entry records what was done, the reasoning behind the decisions, the results, any problems and how they were resolved, and what came next. Commands, tool versions and accessions sufficient for replication are collected in the "Reproducibility" section.*

---

## Phase 1 — Project setup and first reconstruction (7–10 August 2026)

**Aim.** Set the project up and get from a genome to a drawn cassette array.

**What I did.** Created the repository, recorded the IntegronFinder 2.0.2 environment, and wrote the first cassette-array reconstruction from IntegronFinder output, followed by metadata/length annotation and a first visualisation of the array. Recorded notes from the initial supervisor meeting (scope and direction).

**Reasoning.** Pinned the detection environment in the first week so every later result is reproducible, and built the reconstruction step first because every downstream analysis and the visualiser depend on it.

**Next.** Install and compare the detection tools properly.

---

## Phase 2 — Detection tools and tool comparison (27–28 August 2026)

**Aim.** Install the detectors and compare them on a real superintegron.

**What I did.** Installed and recorded HattCI; documented setting I-VIP aside (class-1 scope; no integrons on these genomes). Ran IntegronFinder and HattCI on the *V. cholerae* N16961 chromosome II superintegron, added a HattCI Viterbi-score threshold analysis, wrote a tool-comparison script and statistics, and extracted the cassette table and sequences (59 cassettes from N16961 chr II). Added a zoomed view for inspecting individual cassettes and their attC sites.

**Reasoning.** Compared detectors on a common basis (a shared genome, with a swept HattCI score threshold) rather than at a single arbitrary cut-off, so the comparison would be fair. Recorded the I-VIP decision with its reason rather than silently dropping it.

**Next.** Annotate the cassettes and run the pipeline across more genomes.

---

## Phase 3 — Annotation and multi-genome pipeline (29 August 2026)

**Aim.** Give cassettes real names and automate the whole path.

**What I did.** Added cassette annotation by BLASTx against Swiss-Prot plus CARD/RGI for resistance genes, wrapped detection, reconstruction and annotation into one automated pipeline script, and ran it on four further *Vibrio* genomes (*V. vulnificus* 0, *V. parahaemolyticus* 44, *V. cholerae* O395 50, *A. fischeri* 3 cassettes). Found and removed an earlier result set that had used the wrong genome.

**Reasoning.** The low and highly variable cassette counts across these genomes were the first clear sign that detection settings strongly affect the result — the observation that later motivated the formal benchmark. Corrected the wrong-genome input as soon as it was spotted, to stop it propagating.

**Next.** Scale up and quantify.

---

## Phase 4 — Genus-wide survey and quantitative analyses (31 August 2026)

**Aim.** Put the single-genome observations on a population footing.

**What I did.** Ran the detection pipeline across a genus-wide survey of 988 complete Vibrionaceae genomes. Measured HattCI concordance with IntegronFinder (84%), ran DefenseFinder (3,741 anti-phage systems across 78% of genomes), and carried out a Cayley/permutation analysis of cassette order, which showed order is conserved above ~99% ANI and random below. Recorded all software and database versions and tidied the repository (ignore SLURM logs and regenerable intermediates).

**Reasoning.** The Cayley/permutation result distinguishes cassette-order rearrangement from overall genome relatedness, which whole-genome ANI alone cannot do — this is the quantitative basis for the later distance-network view. Recording versions here fixed the reproducibility baseline for the whole project.

**Gap identified.** HattCI–IntegronFinder concordance varied between genomes, raising the question the final benchmark answered: is the disagreement a threshold effect or a model-level difference?

---

## Phase 5 — Interactive visualiser and distance network (11–18 September 2026)

**Aim.** Turn the static plots into an interactive, attC-resolved viewer with a comparison mode.

**What I did.** Built the v0.1 interactive visualiser (a `.integrons` parser and `plot_integron()`), rendered it on N16961, then a v0.2 multi-track view with gene arrows, attC sites and lengths, and per-cassette function categories. Added multi-genome comparison data, a Cayley/move-to-front distance-network view, and a v2 comparison in which ribbons are shaded by per-pair amino-acid identity. Bundled example `.integrons` and annotation tables as reproducible test cases.

**Consultation.** Reviewed the deliverable shape with J. Engelstädter. Outcome: an R/Shiny app *is* the webpage tool, so "keep it Shiny" and "a webpage where users upload a sequence" are the same deliverable; agreed the input is a genome (the tool runs IntegronFinder itself), attC belongs to its cassette and is revealed on zoom, comparison is pairwise, and annotation comes from public databases.

**Reasoning.** Made the ribbons a similarity gradient instead of a flat shared-family colour so each link shows how conserved it is, and kept ribbons neutral until the identity step runs rather than showing fabricated values.

**Next.** Package the viewer and engine into a usable tool.

---

## Phase 6 — intervis packaged and released (28 September – 9 October 2026)

**Aim.** One engine behind an R/Shiny interface and a CLI, released so others can run it.

**What I did.** Built the engine (run IntegronFinder `--local-max`, reconstruct cassette intervals between consecutive attC sites, name cassette proteins by DIAMOND blastp against Swiss-Prot, emit the viewer data) with the R/Shiny interface and CLI on top. Split the tool into its own repository (`NqynX/intervis`) and made the v1.0.0 public release, then fixed it so a fresh clone runs out of the box — a supervisor clone had shown only the demo genomes and nothing on Run, because the reference genomes were git-ignored and the environment file was incomplete; corrected by committing the reference genomes, completing `environment.yml`, and fixing the backend path resolution.

**Reasoning.** Kept detection as an established method (IntegronFinder) rather than writing a new attC classifier; the contribution is the integron-native end-to-end path with attC-resolved, zoomable cassettes and built-in pairwise comparison.

**Next.** Evaluate the detection step formally against genomes with independent experimental evidence.

---

## Phase 7 — Final benchmark on characterised genomes (early–mid October 2026)

**Aim.** Replace the preliminary comparison with a benchmark on genomes whose integrons are experimentally described, and check the metric is appropriate.

**What I did.** Selected five genomes with experimental evidence (*V. cholerae* 569B, *A. fischeri* DSM 507, *V. vulnificus* ATCC 27562, *V. parahaemolyticus* ATCC 17802, *X. campestris* ATCC 33913). Ran IntegronFinder (default and `--local-max`), HattCI, I-VIP and BacAnt; recorded attC counts, spans and database matches. Scored precision, recall and F1 against the IntegronFinder `--local-max` operational reference (midpoints within 100 bp). Then, questioning whether a single-threshold F1 suits a continuous-score detector, added average precision (AP) — a threshold-independent summary from the precision–recall curve.

**Results and reasoning.** Local search recovered far more attC in the Vibrionaceae arrays than the default (e.g. 45→118 in 569B, 22→126 in *A. fischeri*) but not in *X. campestris*. HattCI agreement was strongly genome-dependent (F1 0.67 to 0.01; AP 0.95 to 0.004). On *A. fischeri*, precision stayed below 0.05 at every threshold, so no cut-off recovered that array — a model-level failure. This matches the model types: IntegronFinder scores attC secondary structure (covariance model), HattCI the primary sequence (HMM), and HattCI fails on the divergent *A. fischeri* repeats described by Rowe-Magnus et al. (2003).

**Inquiry.** F1 measures one operating point and can hide that a tool would be accurate at a different threshold; AP is the better-suited metric, and both are reported.

**Problems and fixes.** BacAnt needed a dependency chain resolved (Biopython 1.76 for a removed import, a BCBio GFF parser patch, several pip/conda installs) before it ran; HattCI runs as a local compiled binary, not an environment command — both recorded so the runs reproduce.

**Next.** Measure genuine accuracy, not concordance.

---

## Phase 8 — Synthetic benchmark (mid October 2026)

**Aim.** Measure absolute recall/precision against a known attC catalogue and test the structure-versus-sequence mechanism directly.

**What I did.** Built synthetic genomes with a known truth: excised the native superintegron from *V. cholerae* 569B (the bare backbone returned no calls, confirming it was clean), reinserted an 80-cassette array rebuilt from the genome's own integrase and real cassette units, and diverged the attC from 0–30% in 5% steps two ways — random substitution (breaks the fold) and compensatory substitution (changes sequence by the same amount while keeping the inverted repeat paired). Ran the four detectors and scored against the known coordinates (one genome per condition, seed 20260101).

**Results.** IntegronFinder `--local-max` recall fell from 1.00 to 0.06 by 15% random divergence, but was 0.48 at 15% and 0.36 at 25% under compensatory divergence; HattCI declined almost identically under both schemes. Structure preservation rescues the covariance-model tool but not the sequence HMM — the controlled version of the *A. fischeri* result.

**Honest limits.** The inserted attC are real 569B sites IntegronFinder first called (a selection bias), so the comparison rests on the slope of recall against divergence and the gap between schemes, not the starting value; complementarity at the chosen positions is a construction, not a validated fold; with one realisation per condition the trends are descriptive, not statistical.

---

## Phase 9 — Report, reproducibility and records (October 2026)

**Aim.** Write up the work within the word limit, make it reproducible, and keep the record.

**What I did.** Drafted and revised the IMRaD report across successive versions, integrating the detection benchmark, the F1/AP concordance analysis and the synthetic benchmark; built the figures and tables; kept the main text under 2,500 words by moving reproducibility detail to an appendix. Addressed supervisor author-queries and a round of editorial-review comments (scoring and AP definitions, single-realisation disclosure, BacAnt output reconciliation, figure-label neutrality). Verified the intervis annotation parameters against the source code before stating them. Reproduced the case-study count (171 attC → 170 intervals for N16961). Published the analysis scripts, figures, benchmark workbook and this record to the repository.

**Consultation.** Incorporated supervisor feedback (J. Engelstädter, A. Letten) across review cycles.

**Future work.** A curated attC reference across integron types for genuine sensitivity estimates; replicated synthetic mutation runs with structural validation; family-level cassette matching; and linking observed arrays to models of integron evolution.

---

## Reproducibility

**Environment.** UQ Bunya HPC, micromamba environments (`integronfinder`, `hattci`, `bacant`, `ivip`); SLURM (account a_uqccr).

**Detection tools.** IntegronFinder 2.0.2 (Prodigal 2.6.3, HMMER 3.3.2, INFERNAL 1.1.4); HattCI 1.0b; I-VIP v1.2; BacAnt 3.4.0.
- IntegronFinder: `integron_finder --circ [--local-max] <genome.fna>`
- HattCI: `./HattCI/hattci.out -b <genome.fna> <out.txt>` (candidates at Viterbi ≥ 2)
- I-VIP, BacAnt: bundled databases, default settings

**intervis annotation.** DIAMOND 2.2.5 blastp of translated cassette proteins against UniProtKB/Swiss-Prot (575,748 sequences, downloaded 11 June 2026); naming step `diamond blastp --more-sensitive --max-target-seqs 1 --evalue 1e-5` (top hit per protein); pairwise comparison `diamond blastp --more-sensitive -k 50 --evalue 1e-5`, reciprocal best hits at ≥30% identity and ≥0.5 coverage (alignment length / shorter protein).

**Scoring.** A call matches a reference site when midpoints are within 100 bp on the same replicon; matching is one-to-one (greedy nearest, each reference used once). Precision = matched/all calls; recall = matched/all reference sites; F1 = harmonic mean. AP = mean precision recorded at each recovered site when calls are ranked by score (a step summary of the precision–recall curve, not the trapezoidal area).

**Synthetic generator.** `make_synthetic.py` (seed 20260101, 80 attC per array, 0–30% divergence in 5% steps, naive and compensatory schemes, 7 conserved bases fixed at each attC end); `run_synth_tools.sh`; `score_synth.py`.

**Accessions.** Reference assemblies in Table 1 of the report; N16961 visualisation GCF_000006745.1 (chromosome II NC_002506.1).

## Version-control history

Two repositories hold the project. The **vibrio-integron-pipeline** development repo carries the day-to-day work: 40 commits from the initial commit on 7 Aug 2026 through the annotation tables on 18 Sep 2026 (project setup and reconstruction; HattCI install and tool comparison; annotation and the automated pipeline; the 988-genome genus survey, concordance and Cayley analysis; and the interactive visualiser and distance network). The tool was then split into the **NqynX/intervis** release repo: v1.0.0 on 28 Sep 2026, the out-of-the-box fix on 9 Oct 2026 (`6cb0bca`, `git describe` v1.0.0-3-g6cb0bca), and the analysis, figures and records added on 14 Oct 2026. Per-commit notes are recorded against each commit in both repositories.

## Supporting materials (graphics and data)

- **Table 1** — reference genomes and experimental evidence.
- **Table 2** — detection counts across the five genomes.
- **Table 3** — site-level concordance (F1 and AP).
- **Figure 1** — integron architecture schematic.
- **Figure 2** — attC detected per tool, with reported-cassette reference lines.
- **Figure 3** — synthetic recall vs attC divergence (random vs compensatory).
- **Figure A1** — HattCI precision–recall curves (average precision).
- **Spreadsheet** — `SCIE3220_detection_benchmark.xlsx` (detection, F1, AP, synthetic sheets).
- **Scripts** — `score_ap.py`, `make_synthetic.py`, `run_synth_tools.sh`, `score_synth.py`; figure scripts; intervis engine/CLI/app in the repository.

*Note for submission: add the exact dates of supervisor meetings where marked, and confirm the phase dates against your own calendar.*
