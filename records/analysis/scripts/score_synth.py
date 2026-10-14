#!/usr/bin/env python3
# =============================================================================
# score_synth.py  -  score detectors against the KNOWN synthetic attC truth
# Run from tool-benchmark/ after run_synth_tools.sh. Reads synthetic/*.truth.tsv
# and synth_runs/. For each genome and tool it reports recall, precision and FDR
# against the true inserted attC coordinates (absolute accuracy, not concordance).
# Pure standard library.  Usage:  python3 score_synth.py
# =============================================================================
import glob, os

SYN = "synthetic"
OUT = "synth_runs"
TOL = 100          # bp: a call matches a true attC if midpoints are within TOL
VITERBI = 2.0      # HattCI score threshold (same as the main benchmark)

# --------------------------------------------------------------- call parsers
def if_mids(outdir):
    out = []
    for p in glob.glob(os.path.join(outdir, "Results_Integron_Finder_*", "*.integrons")):
        for ln in open(p):
            if ln.startswith("#") or ln.lower().startswith("id_integron"): continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 8 or f[7] != "attC": continue
            try: out.append((int(f[3]) + int(f[4])) / 2.0)
            except ValueError: pass
    return out

def hattci_mids(path, vmin=VITERBI):
    out = []
    if not os.path.exists(path): return out
    for ln in open(path):
        f = ln.split()
        if len(f) < 5 or f[0] == "hit": continue
        try:
            s, e, v = int(f[2]), int(f[3]), float(f[4])
        except ValueError:
            continue
        if v >= vmin: out.append((s + e) / 2.0)
    return out

def bacant_mids(outdir):
    # BacAnt is class-1 oriented; parse any integron coordinate table if present.
    out = []
    for p in glob.glob(os.path.join(outdir, "**", "*.tsv"), recursive=True):
        base = os.path.basename(p).lower()
        if "integron" not in base: continue
        with open(p) as fh:
            head = fh.readline().lower().split("\t")
            def col(*names):
                for n in names:
                    if n in head: return head.index(n)
                return None
            cs, ce = col("start", "attc_start", "pos_beg"), col("end", "attc_end", "pos_end")
            if cs is None or ce is None: continue
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                try: out.append((int(f[cs]) + int(f[ce])) / 2.0)
                except (ValueError, IndexError): pass
    return out

# --------------------------------------------------------------- matching/score
def score(truth, calls):
    t = sorted(truth); used = [False] * len(t); tp = 0
    for c in sorted(calls):
        best, bd = -1, TOL + 1
        for i, m in enumerate(t):
            if used[i]: continue
            d = abs(m - c)
            if d <= TOL and d < bd: bd, best = d, i
        if best >= 0: used[best] = True; tp += 1
    ncall = len(calls); ntrue = len(t)
    recall = tp / ntrue if ntrue else float("nan")
    prec = tp / ncall if ncall else float("nan")
    fdr = 1 - prec if ncall else float("nan")
    return ntrue, ncall, tp, recall, prec, fdr

def load_truth(g):
    p = os.path.join(SYN, g + ".truth.tsv")
    return [float(l.split("\t")[1]) for l in open(p).read().splitlines()[1:]] if os.path.exists(p) else []

# --------------------------------------------------------------- run
man = os.path.join(SYN, "manifest.tsv")
rows = [l.split("\t") for l in open(man).read().splitlines()[1:]] if os.path.exists(man) else []

TOOLS = ["IFdef", "IFloc", "HattCI", "BacAnt"]
def tool_calls(tool, g):
    if tool == "IFdef":  return if_mids(os.path.join(OUT, "IFdef_" + g))
    if tool == "IFloc":  return if_mids(os.path.join(OUT, "IFloc_" + g))
    if tool == "HattCI": return hattci_mids(os.path.join(OUT, "hattci_%s.txt" % g))
    if tool == "BacAnt": return bacant_mids(os.path.join(OUT, "bacant_" + g))
    return []

print("# per-genome scores vs known truth (TOL=%d bp, HattCI Viterbi>=%g)" % (TOL, VITERBI))
print("%-26s %-8s %5s %6s %5s %6s %6s %6s" %
      ("genome", "tool", "nTrue", "nCall", "TP", "Recall", "Prec", "FDR"))
print("-" * 80)
agg = {}   # (mode, lvl, tool) -> (recall, prec)
backbone_fp = {}
for r in rows:
    g, lvl, mode = r[0], int(r[1]), r[2]
    truth = load_truth(g)
    for tool in TOOLS:
        calls = tool_calls(tool, g)
        ntrue, ncall, tp, rec, prec, fdr = score(truth, calls)
        if mode == "backbone":
            backbone_fp[tool] = ncall        # all calls here are false positives
            print("%-26s %-8s %5d %6d %5d %6s %6s %6s" %
                  (g, tool, ntrue, ncall, tp, "-", "-", "-"))
        else:
            agg[(mode, lvl, tool)] = (rec, prec)
            print("%-26s %-8s %5d %6d %5d %6.3f %6.3f %6.3f" %
                  (g, tool, ntrue, ncall, tp, rec, prec, fdr))

print("\n# backbone negative control - false-positive calls in integron-free DNA")
for tool in TOOLS:
    print("  %-8s %d calls" % (tool, backbone_fp.get(tool, 0)))

# ---- paste-back block: recall vs divergence, per mode and tool ----
levels = sorted({int(r[1]) for r in rows if r[2] != "backbone"})
print("\n# ==== RECALL vs DIVERGENCE (paste this block back) ====")
print("mode,tool,divergence,recall,precision")
for mode in ("naive", "structure"):
    for tool in TOOLS:
        for lvl in levels:
            if (mode, lvl, tool) in agg:
                rec, prec = agg[(mode, lvl, tool)]
                print("%s,%s,%d,%.3f,%.3f" % (mode, tool, lvl, rec, prec))
