# Precision-Recall curve + Average Precision (AP / PR-AUC) for the attC detectors,
# scored against the IntegronFinder --local-max reference (attC midpoints <=100 bp).
# Threshold-independent: each tool's own score is swept (HattCI = Viterbi; IF default = E-value).
# Runs on the files already in /scratch/user/uqcngu19/tool-benchmark/ ; stdlib only.
import glob, os
GENOMES = ["vc569B","afDSM","vv275","vp178","xcc"]
TOL = 100

def ref_mids(g):
    out = {}
    p = glob.glob("runs/IF_local_%s/Results_Integron_Finder_%s/%s.integrons" % (g,g,g))
    if not p: return out
    for ln in open(p[0]):
        if ln.startswith("#") or ln.startswith("ID_integron"): continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 8 or f[7] != "attC": continue
        out.setdefault(f[1], []).append((int(f[3])+int(f[4]))/2.0)
    return out

def hattci_scored(g):
    # returns list of (replicon, midpoint, score) ; higher score = more confident
    out = []
    p = "runs/hattci_%s.txt" % g
    if not os.path.exists(p): return out
    for ln in open(p):
        f = ln.split()
        if len(f) < 5 or f[0] == "hit": continue
        try: s,e,v = int(f[2]),int(f[3]),float(f[4])
        except: continue
        out.append((f[1], (s+e)/2.0, v))
    return out

def ifdefault_scored(g):
    # IntegronFinder default calls, ranked by E-value (smaller = more confident -> use -log sort)
    out = []
    p = glob.glob("runs/IF_nolocal_%s/Results_Integron_Finder_%s/%s.integrons" % (g,g,g))
    if not p: return out
    for ln in open(p[0]):
        if ln.startswith("#") or ln.startswith("ID_integron"): continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 8 or f[7] != "attC": continue
        try: ev = float(f[6])
        except: ev = 1.0
        score = -ev                      # higher = better (smaller E-value)
        out.append((f[1], (int(f[3])+int(f[4]))/2.0, score))
    return out

def average_precision(ref, cand):
    """cand = list of (replicon, mid, score). Returns AP, curve=[(recall,precision)], nref, maxrecall."""
    nref = sum(len(v) for v in ref.values())
    if nref == 0: return 0.0, [], 0, 0.0
    # mutable matched flags per replicon
    refs = {k: [[m, False] for m in sorted(v)] for k, v in ref.items()}
    cand = sorted(cand, key=lambda x: x[2], reverse=True)   # high score first
    tp = 0; fp = 0; ap = 0.0; curve = []
    for rep, mid, sc in cand:
        rs = refs.get(rep, [])
        best = -1; bestd = TOL + 1
        for i, (rm, used) in enumerate(rs):
            if used: continue
            d = abs(rm - mid)
            if d <= TOL and d < bestd: bestd = d; best = i
        if best >= 0:
            rs[best][1] = True; tp += 1
            prec = tp / (tp + fp)
            ap += prec                      # sum precision at each true-positive event
            curve.append((tp / nref, prec))
        else:
            fp += 1
    ap = ap / nref
    maxrec = tp / nref
    return ap, curve, nref, maxrec

print("%-8s %6s %7s %7s %9s %9s" % ("genome","nRef","nHat","HatAP","HatMaxRec","IFdefAP"))
print("-"*52)
allcurves = {}
for g in GENOMES:
    ref = ref_mids(g)
    hc = hattci_scored(g)
    ap_h, curve_h, nref, maxrec = average_precision(ref, hc)
    ap_i, _, _, _ = average_precision(ref, ifdefault_scored(g))
    allcurves[g] = curve_h
    print("%-8s %6d %7d %7.3f %9.3f %9.3f" % (g, nref, len(hc), ap_h, maxrec, ap_i))

# write full curves + print a downsampled version for plotting
with open("pr_curve_hattci.csv","w") as fh:
    fh.write("genome,recall,precision\n")
    for g,c in allcurves.items():
        for r,p in c: fh.write("%s,%.5f,%.5f\n" % (g,r,p))
print("\n# downsampled HattCI PR points (paste this block back):")
print("genome,recall,precision")
for g,c in allcurves.items():
    if not c: 
        print("%s,0,0" % g); continue
    n = len(c); step = max(1, n//18)
    pts = c[::step]
    if pts[-1] != c[-1]: pts.append(c[-1])
    for r,p in pts: print("%s,%.4f,%.4f" % (g,r,p))
