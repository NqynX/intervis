# SCIE3220 Research Project — Record Keeping

**Student:** Liam Nguyen (48596804)
**Project:** Development of a bioinformatics pipeline for identifying and visualising integron cassette arrays in *Vibrio* genomes (intervis)
**Supervisors:** Jan Engelstädter; Andrew Letten (Engelstädter Lab)
**Code and data:** https://github.com/NqynX/intervis
**Compute:** UQ Bunya HPC cluster (micromamba environments, SLURM)

*This is a chronological record of a computational/bioinformatics project. Entries are grouped by work phase. Dates are to the week; where an exact date matters it is given. Each entry records what was done, the reasoning behind the decisions, the results, any problems and how they were resolved, and what came next. Commands, tool versions and accessions sufficient for replication are collected in the "Reproducibility" section at the end.*

---

## Phase A — Background, scope and aims (early September 2026)

**Aim.** Define a project around integrons in the Vibrionaceae and settle a concrete, assessable deliverable.

**What I did.** Read the integron literature (Hall and Collis 1995; Mazel et al. 1998; Rowe-Magnus et al. 2001, 2003; Cury et al. 2016) and drafted the scientific introduction. Integrons capture and rearrange gene cassettes by site-specific recombination; in the Vibrionaceae they form chromosomal superintegrons that can exceed 100 kb and hold more than a hundred cassettes. The cassette boundaries are the attC recombination sites, which are conserved mainly through structure and vary widely in sequence, so detecting them is the technical bottleneck for reconstructing an array.

**Decisions and reasoning.** Set two aims: (1) test how the choice of detection tool and its settings changes the reconstructed array, and (2) build intervis, a reproducible pipeline that turns a detection result into an annotated, interactive, attC-resolved view with a built-in pairwise comparison. The first aim is needed to justify the detection step the second aim depends on — I did not want to build a visualiser on top of a detection step I had not evaluated.

**Next.** Build the viewers and confirm the detection tooling.

---

## Phase B — Prototype viewers and cassette-order distance (September 2026)

**Aim.** Prove the interactive visualisation path on real *Vibrio* data and prototype a way to compare cassette order between arrays.

**What I did.** Built two browser viewers in HTML/JS on real *Vibrio* superintegron data: a single-array view ("Cartographer") and a stacked comparison view ("Synteny") that links cassettes of the same protein family between genomes with ribbons. Added a distance-network prototype: each superintegron is reduced to an ordered list of cassette-family identifiers, and pairs are compared on the shared single-copy backbone using Cayley (transposition) distance and a move-to-front distance, with a random-permutation baseline and a Mash whole-genome distance carried alongside for contrast.

**Decisions and reasoning.** Made the comparison ribbons a similarity gradient (per-pair amino-acid identity) instead of a flat "shared family" colour, because a flat colour hides how conserved each link actually is. Kept the ribbons neutral until the identity step is run rather than showing fabricated values. Chose Cayley/move-to-front distance because they measure rearrangement of the cassette order specifically, which Mash/ANI does not — ANI answers "how related are the genomes", the permutation distance answers "how rearranged is the array".

**Consultation.** Reviewed the deliverable shape with J. Engelstädter. Outcome: an R/Shiny app *is* the webpage tool, so "keep it Shiny" and "a webpage where users upload a sequence" are the same thing. Agreed the input should be a genome sequence (intervis runs IntegronFinder itself), attC should belong to its cassette and be revealed on zoom, comparison is pairwise in the web tool, and annotation should come from public databases.

**Next.** Fold IntegronFinder and annotation into a single engine behind the viewers.

---

## Phase C — intervis engine and Shiny interface (late September – October 2026)

**Aim.** One engine, two front-ends: a Python engine shared by an R/Shiny web interface and a command-line tool.

**What I did.** Built the engine: run IntegronFinder with `--local-max`, parse the `.integrons` output by column name, orient the array with the integrase on the left, reconstruct cassette intervals between consecutive attC sites, translate the cassette proteins and name them by DIAMOND blastp against Swiss-Prot, and emit the data the viewers draw from. Added the pairwise comparison: all cassette proteins searched against one another with DIAMOND, reciprocal best hits between arrays retained above identity and coverage thresholds, each link shaded by its own identity. Wrapped the engine in an R/Shiny interface that uploads one or two genomes and embeds the viewers. Iterated the GUI over many local snapshots.

**Decisions and reasoning.** Kept detection as an established method (IntegronFinder) rather than writing a new attC classifier — the contribution is the integron-native end-to-end path, not a new detector. Defined cassette intervals only between consecutive attC sites, so an array with *n* attC gives *n*−1 intervals; the promoter-proximal region between attI and the first attC is not counted as a cassette interval, which I noted so the interval count is not mistaken for a biological cassette total. Labelled proteins with no database hit as hypothetical rather than guessing, since a high hypothetical fraction is a real property of superintegron cassettes.

**Problems and fixes.** Running IntegronFinder from the tool needs a host with the `integronfinder` environment installed; resolved by running on Bunya in a micromamba environment. A public deployment would need a server with that environment and job handling — parked as future work, not a Phase-1 concern.

**Next.** Evaluate the detection step on genomes with independent experimental evidence.

---

## Phase D — Reference genomes and validation framework (early–mid October 2026)

**Aim.** Choose genomes whose integrons are experimentally described, so detection can be judged against something other than another tool.

**What I did.** Selected five genomes: four Vibrionaceae superintegrons (*V. cholerae* 569B, *A. fischeri* DSM 507, *V. vulnificus* ATCC 27562, *V. parahaemolyticus* ATCC 17802) and one shorter *Xanthomonas campestris* ATCC 33913 integron with unusually short attC sites. Recorded accessions and the specific experimental evidence for each (Table 1).

**Decisions and reasoning.** No complete, base-resolved experimental catalogue of attC sites exists for these assemblies, so true sensitivity cannot be measured directly. Chose the IntegronFinder `--local-max` call set as an *operational* reference because it gave the most complete reconstructions, and was explicit throughout that scores against it measure concordance, not absolute accuracy. Included *X. campestris* as a structurally contrasting, short-attC integron so the benchmark was not only long Vibrionaceae arrays.

**Gap identified.** The operational reference is itself a tool, so concordance is not accuracy — this gap is what later motivated the synthetic benchmark (Phase G).

---

## Phase E — Detection benchmark (mid October 2026)

**Aim.** Measure how detector choice and search depth change the reconstructed array across the five genomes.

**What I did.** Ran IntegronFinder 2.0.2 twice per genome (default and `--local-max`), HattCI 1.0b on both strands (Viterbi ≥ 2), and I-VIP v1.2 and BacAnt 3.4.0 with their bundled databases. Recorded attC counts, the integrase-anchored span, and database matches (Table 2; Figure 2).

**Results.** Local search increased attC counts in all four Vibrionaceae genomes (e.g. 45→118 in 569B, 22→126 in *A. fischeri*) but not in *X. campestris* (22 either way). The `--local-max` *V. vulnificus* array spanned 156.8 kb, close to the 158.3-kb superintegron resolved in the complete genome of the same strain. I-VIP returned no integrons; BacAnt returned only 0–4 class-1 database matches.

**Problems and fixes.** BacAnt would not run out of the box: resolved a dependency chain (matplotlib, reportlab, colorlog, networkx; `dna_features_viewer` installed with `pip --no-deps`), downgraded Biopython to 1.76 for a removed `Bio.Alphabet` import, and patched a `SeqFeature.SimpleLocation` call in the BCBio GFF parser. After that BacAnt ran cleanly. HattCI runs as a local compiled binary (`./HattCI/hattci.out`), not an environment command — noted this so the runs are reproducible.

**Reasoning.** Because each cassette is bounded by attC sites, a missed site merges two cassettes, so search depth changes the inferred architecture, not only the site count — this is why the search setting matters for the downstream comparison.

**Next.** Quantify how closely the other tools agree with the operational reference.

---

## Phase F — Concordance and threshold-independent accuracy (mid–late October 2026)

**Aim.** Score agreement with the operational reference, and check the metric is appropriate.

**What I did.** Scored precision, recall and F1 with `score_f1.py` (a call matches a reference site when midpoints are within 100 bp, one-to-one on the same replicon) for the detectors with comparable attC coordinates (Table 3). Then questioned whether a single-threshold F1 is the right summary for HattCI, which produces a continuous score. Added average precision (AP) with `score_ap.py`: rank HattCI candidates by Viterbi score, sweep every cut-off, and take the mean precision at each recovered site — a threshold-independent summary (Figure A1).

**Results and reasoning.** HattCI agreement was strongly genome-dependent (F1 0.67 in 569B to 0.01 in *A. fischeri*; AP 0.95 in *X. campestris* to 0.004 in *A. fischeri*). On *X. campestris* HattCI held precision 1.00 up to ~95% recall, so a stricter threshold would have been nearly exact there; on *A. fischeri* precision stayed below 0.05 at every threshold, so no cut-off recovered that array. This confirmed the model-level reading: IntegronFinder uses a covariance model that scores attC secondary structure, HattCI an HMM that scores the primary sequence, and HattCI fails on the divergent *A. fischeri* repeats because its model is sequence-based (consistent with the species-specific repeat divergence reported by Rowe-Magnus et al. 2003).

**Inquiry.** Recognised that F1 measures one operating point and can hide that a tool would be accurate at a different threshold; AP is the better-suited, threshold-independent metric, and the two are reported together. Confirmed HattCI's model type from its publication (hidden Markov model) rather than assuming.

**Next.** Replace concordance with a genuine accuracy measurement.

---

## Phase G — Synthetic benchmark with known truth (late October 2026)

**Aim.** Measure absolute recall/precision against a known attC catalogue, and test the structure-versus-sequence mechanism directly.

**What I did.** Built synthetic genomes with `make_synthetic.py`: excised the native superintegron from the *V. cholerae* 569B chromosome (IntegronFinder `--local-max` span plus a 5-kb margin; the bare backbone returned no attC calls, confirming it was clean), and reinserted an 80-cassette array rebuilt from the genome's own integrase and real cassette units, recording every attC coordinate. Diverged the attC sequences from 0 to 30% in 5% steps two ways: random substitution (breaks the fold) and compensatory substitution on the two arms (changes sequence by the same amount while keeping the inverted repeat paired). Ran the four detectors (`run_synth_tools.sh`) and scored against the known coordinates (`score_synth.py`). One genome per condition, random seed 20260101.

**Results.** On the empty backbone, IntegronFinder and BacAnt made no calls and HattCI eight. With the array present, IntegronFinder `--local-max` recall fell from 1.00 to 0.06 by 15% *random* divergence, but was 0.48 at 15% and 0.36 at 25% under *compensatory* divergence of the same magnitude; HattCI declined almost identically under both schemes (0.20 and 0.25 at 15%). So structure preservation rescues the covariance-model tool but not the sequence HMM — the controlled version of the *A. fischeri* result (Figure 3).

**Inquiry and honest limits.** Designed this specifically to separate sequence change from structure loss, which the real genomes cannot do. Recorded the caveats plainly: the inserted attC are real 569B sites that IntegronFinder first called (a selection bias), so the comparison rests on the *slope* of recall against divergence and the gap between the two schemes, not the starting value; complementarity at the chosen positions is a construction, not an experimentally validated fold; and with one realisation per condition the trends are descriptive, not statistical. Replicated mutation runs with structural validation would test the mechanism more directly.

**Next.** Write up the benchmark and the tool, and publish the code.

---

## Phase H — Report writing and revision (2–13 October 2026)

**Aim.** An IMRaD report within the word limit, accurate to the results, readable, and defensible.

**What I did.** Drafted and revised the report through successive versions. Integrated the detection benchmark, the F1/AP concordance analysis and the synthetic benchmark; built the figures (detection counts, precision–recall/AP, synthetic recall) and the tables; kept the main text under the 2500-word limit by moving reproducibility detail to an appendix. Addressed supervisor author-queries (tool versions, database release, DIAMOND and classification settings) and a round of editorial-review comments (scoring definition, average-precision definition as a step summary rather than a trapezoidal area, single-realisation disclosure, BacAnt output reconciliation, figure-label neutrality, rounding of the *A. fischeri* AP).

**Decisions and reasoning.** Kept the claims hedged to what the evidence supports: concordance is not accuracy; a count close to the published cassette number is completeness, not correctness (shown for *V. vulnificus* and 569B, with no reference line drawn where only partial cassette sets have been cloned); and detector differences are explained mechanistically but not over-claimed as the proven cause of the *A. fischeri* result. Verified the intervis annotation parameters against the source code rather than from memory before stating them, which corrected a draft that had merged the Swiss-Prot naming step with the comparison step.

**Consultation.** Incorporated supervisor feedback across review cycles (J. Engelstädter, A. Letten).

---

## Phase I — Publishing, reproducibility and verification (October 2026)

**Aim.** Make the work reproducible by a reviewer from the repository alone.

**What I did.** Published intervis to GitHub (NqynX/intervis), with the engine, CLI, R/Shiny app, environment file and bundled reference *Vibrio* genomes so a supervisor can clone and run the stored examples. Confirmed the release at commit 6cb0bca (`git describe` v1.0.0-3-g6cb0bca). Recorded the exact tool versions, commands, database release and accessions for the appendix. Independently reproduced the case-study count: an IntegronFinder `--local-max` run on the *V. cholerae* N16961 assembly (GCF_000006745.1) returned 171 attC sites, giving the 170 reconstructed intervals reported.

**Problems and fixes.** An early clone by the supervisor showed only the demo genomes and nothing on "Run", because the real reference genomes were git-ignored and the environment file was empty; fixed by committing the reference genomes and a complete environment file and re-pushing.

**Next / future work.** A curated attC reference spanning integron types for genuine sensitivity estimates; replicated synthetic mutation runs with structural validation; family-level (not one-to-one) cassette matching; and linking observed arrays to models of integron evolution.

---

## Reproducibility

**Environment.** UQ Bunya HPC, micromamba environments (`integronfinder`, `hattci`, `bacant`, `ivip`); SLURM (account a_uqccr).

**Detection tools.** IntegronFinder 2.0.2 (Prodigal 2.6.3, HMMER 3.3.2, INFERNAL 1.1.4); HattCI 1.0b; I-VIP v1.2; BacAnt 3.4.0.
- IntegronFinder: `integron_finder --circ [--local-max] <genome.fna>`
- HattCI: `./HattCI/hattci.out -b <genome.fna> <out.txt>` (candidates at Viterbi ≥ 2)
- I-VIP, BacAnt: bundled databases, default settings

**intervis annotation.** DIAMOND 2.2.5 blastp of translated cassette proteins against UniProtKB/Swiss-Prot (575,748 sequences, downloaded 11 June 2026); naming step `diamond blastp --more-sensitive --max-target-seqs 1 --evalue 1e-5` (top hit per protein); pairwise comparison `diamond blastp --more-sensitive -k 50 --evalue 1e-5`, reciprocal best hits retained at ≥30% identity and ≥0.5 coverage (alignment length / shorter protein).

**Scoring.** A call matches a reference site when midpoints are within 100 bp on the same replicon; matching is one-to-one (greedy nearest, each reference used once). Precision = matched/all calls; recall = matched/all reference sites; F1 = harmonic mean. AP = mean precision recorded at each recovered site when calls are ranked by score (a step summary of the precision–recall curve, not the trapezoidal area).

**Synthetic generator.** `make_synthetic.py` (seed 20260101, 80 attC per array, 0–30% divergence in 5% steps, naive and compensatory schemes, 7 conserved bases fixed at each attC end); `run_synth_tools.sh`; `score_synth.py`.

**Accessions.** Reference assemblies in Table 1 of the report; N16961 visualisation GCF_000006745.1 (chromosome II NC_002506.1).

**Code.** https://github.com/NqynX/intervis (commit 6cb0bca, tag v1.0.0).

**Version-control history (NqynX/intervis).**
- `2d02ba9` (2026-09-28) — Initial public release: intervis v1.0.0 (engine, CLI, R/Shiny app).
- `239b2f0` (2026-09-28) — docs: point citation and changelog links at the public repository.
- `42e42ab` (2026-09-28) — docs: rewrite README in plain style.
- `6cb0bca` (2026-10-09) — bundle three reference *Vibrio* examples so a fresh clone runs out of the box; add biopython and R/Shiny to the environment; fix backend path resolution. (`git describe` v1.0.0-3-g6cb0bca; used for the reproducibility appendix.)

---

## Supporting materials (graphics and data)

- **Table 1** — reference genomes and experimental evidence.
- **Table 2** — detection counts across the five genomes.
- **Table 3** — site-level concordance (F1 and AP).
- **Figure 1** — integron architecture schematic.
- **Figure 2** — attC detected per tool, with reported-cassette reference lines.
- **Figure 3** — synthetic recall vs attC divergence (random vs compensatory).
- **Figure A1** — HattCI precision–recall curves (average precision).
- **Spreadsheet** — `SCIE3220_detection_benchmark.xlsx` (detection, F1, AP, synthetic sheets).
- **Scripts** — `score_f1.py`, `score_ap.py`, `make_synthetic.py`, `run_synth_tools.sh`, `score_synth.py`; figure scripts; intervis engine/CLI/app in the repository.

*Note for submission: dates above are to the week of work; please align them to your own calendar and add the exact dates of supervisor meetings before submitting.*
