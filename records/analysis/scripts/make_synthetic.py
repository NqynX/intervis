#!/usr/bin/env python3
# =============================================================================
# make_synthetic.py  -  build synthetic integron genomes with known attC truth
# -----------------------------------------------------------------------------
# Idea: take a real Vibrio replicon, cut out its native superintegron to get an
# integron-free backbone, then rebuild an array from that genome's own real
# cassette units (ORF + attC) and insert it back. Because we place every attC
# ourselves, we know the exact number and coordinates -> a real gold standard,
# so recall / precision / FDR become absolute, not agreement with another tool.
#
# The attC repeats are then diverged in sequence at a series of levels, two ways:
#   naive      - random substitutions that ignore the fold
#   structure  - substitutions on one arm mirrored onto the other, so the
#                imperfect inverted repeat (the attC secondary structure) is
#                preserved while the primary sequence diverges
# This titrates the structure-vs-sequence hypothesis directly: a structure-aware
# detector (IntegronFinder covariance model) should tolerate more divergence than
# a primary-sequence HMM (HattCI) before recall collapses.
#
# Dependencies: none (plain Python 3). Run on Bunya where the genome + IF output
# live. Inputs it needs (defaults match the tool-benchmark layout):
#   --fna        donor genome FASTA              (runs input .fna)
#   --integrons  IntegronFinder --local-max .integrons for that genome
#   --out        output directory
# =============================================================================
import argparse, os, random, re

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def rc(s): return s.translate(COMP)[::-1]

# ----------------------------------------------------------------- FASTA I/O
def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name is not None: seqs[name] = "".join(buf)
                name = ln[1:].split()[0]; buf = []
            else:
                buf.append(ln.strip())
    if name is not None: seqs[name] = "".join(buf)
    return seqs

def write_fasta(path, name, seq, width=70):
    with open(path, "w") as fh:
        fh.write(">%s\n" % name)
        for i in range(0, len(seq), width):
            fh.write(seq[i:i+width] + "\n")

# ----------------------------------------------- parse IntegronFinder .integrons
def parse_integrons(path):
    """Return (replicon_id, intI(beg,end) or None, [attC (beg,end) sorted])."""
    attc = []; intI = None; replicon = None
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or ln.lower().startswith("id_integron"): continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9: continue
            rep = f[1]
            try: beg, end = int(f[3]), int(f[4])
            except ValueError: continue
            type_elt = f[7]; annot = f[8]
            if replicon is None: replicon = rep
            if rep != replicon: continue          # keep the first replicon only
            if type_elt == "attC":
                attc.append((min(beg, end), max(beg, end)))
            elif annot == "intI" or type_elt == "intI":
                intI = (min(beg, end), max(beg, end))
    attc.sort()
    return replicon, intI, attc

# ----------------------------------------------------------- divergence engines
def _sub(base, rng):
    alt = [b for b in "ACGT" if b != base.upper()]
    return rng.choice(alt) if alt else base

def diverge_naive(attc, pct, core, rng):
    s = list(attc); L = len(s)
    mutable = [i for i in range(core, L - core)]
    n = int(round(pct / 100.0 * L))
    for i in rng.sample(mutable, min(n, len(mutable))):
        s[i] = _sub(s[i], rng)
    return "".join(s)

def diverge_structure(attc, pct, core, rng):
    """Mutate the left arm, mirror the complement onto the right arm, so the
    stem stays paired (structure preserved) while sequence changes. The central
    unpaired region is mutated freely."""
    s = list(attc); L = len(s)
    lo, hi = core, L - core                 # mutable window
    mid = (lo + hi) // 2
    left = list(range(lo, mid))             # 5' arm positions
    changes_target = int(round(pct / 100.0 * L))
    changed = 0
    rng.shuffle(left)
    for i in left:
        j = L - 1 - i                       # paired partner on the 3' arm
        if j <= i: continue
        nb = _sub(s[i], rng)
        s[i] = nb; s[j] = COMP_CHAR(nb)     # compensatory: keep complementarity
        changed += 2
        if changed >= changes_target: break
    # top up with free changes in the central loop if target not reached
    if changed < changes_target:
        loop = [k for k in range(mid, hi) if (L - 1 - k) < lo or (L-1-k) >= hi]
        rng.shuffle(loop)
        for k in loop:
            s[k] = _sub(s[k], rng); changed += 1
            if changed >= changes_target: break
    return "".join(s)

def COMP_CHAR(b):
    return {"A":"T","T":"A","G":"C","C":"G"}.get(b.upper(), b)

def pct_diff(a, b):
    n = sum(1 for x, y in zip(a, b) if x.upper() != y.upper())
    return 100.0 * n / len(a) if a else 0.0

# --------------------------------------------------------------------- builder
def build(args):
    rng = random.Random(args.seed)
    os.makedirs(args.out, exist_ok=True)
    seqs = read_fasta(args.fna)
    replicon, intI, attc = parse_integrons(args.integrons)
    if replicon not in seqs:
        # fall back: longest contig
        replicon = max(seqs, key=lambda k: len(seqs[k]))
    genome = seqs[replicon]
    if len(attc) < 3:
        raise SystemExit("need >=3 attC sites in the donor .integrons")

    # --- native array extent (what we excise) ---
    first_beg = attc[0][0]; last_end = attc[-1][1]
    head_beg = min(first_beg, intI[0]) if intI else first_beg
    ex_lo = max(0, head_beg - args.margin)
    ex_hi = min(len(genome), last_end + args.margin)
    backbone = genome[:ex_lo] + genome[ex_hi:]
    insert_at = ex_lo                                   # reinsert at the cut

    # --- head (integrase + attI region) taken from the real genome ---
    head = genome[intI[0]-1:first_beg-1] if intI else ""

    # --- cassette units: start of attC_i .. start of attC_{i+1} ---
    units = []           # (unit_seq, attc_len)  attC sits at the 5' end of unit
    for i in range(len(attc) - 1):
        u_beg = attc[i][0]; u_end = attc[i+1][0]        # 0-based half-open below
        unit = genome[u_beg-1:u_end-1]
        aL = attc[i][1] - attc[i][0] + 1
        if len(unit) > aL + 10:                         # need attC + some spacer
            units.append((unit, aL))
    if not units:
        raise SystemExit("could not extract cassette units")
    K = min(args.cassettes, len(units))
    units = units[:K]

    levels = [int(x) for x in args.levels.split(",")]
    manifest = []

    # negative control: bare backbone, no array
    bb_name = "synthetic_backbone"
    write_fasta(os.path.join(args.out, bb_name + ".fna"), bb_name, backbone)
    open(os.path.join(args.out, bb_name + ".truth.tsv"), "w").write("replicon\tattc_mid\n")
    manifest.append((bb_name, 0, "backbone", 0))

    for mode in ("naive", "structure"):
        divg = diverge_naive if mode == "naive" else diverge_structure
        for lvl in levels:
            arr = []; truth = []; obs = []
            for (unit, aL) in units:
                attc_seq = unit[:aL]; spacer = unit[aL:]
                mut = attc_seq if lvl == 0 else divg(attc_seq, lvl, args.core, rng)
                obs.append(pct_diff(attc_seq, mut))
                # position of this attC midpoint within the assembled array
                mid_in_arr = len(head) + sum(len(x) for x in arr) + aL / 2.0
                truth.append(mid_in_arr)
                arr.append(mut + spacer)
            array_seq = head + "".join(arr)
            synth = backbone[:insert_at] + array_seq + backbone[insert_at:]
            # truth coords shift by insert_at
            name = "synthetic_%s_d%02d" % (mode, lvl)
            write_fasta(os.path.join(args.out, name + ".fna"), name, synth)
            with open(os.path.join(args.out, name + ".truth.tsv"), "w") as fh:
                fh.write("replicon\tattc_mid\n")
                for m in truth:
                    fh.write("%s\t%d\n" % (name, int(insert_at + m)))
            meanobs = sum(obs) / len(obs) if obs else 0.0
            manifest.append((name, lvl, mode, round(meanobs, 2)))

    with open(os.path.join(args.out, "manifest.tsv"), "w") as fh:
        fh.write("genome\ttarget_divergence\tmode\tmean_observed_divergence\tn_attC\n")
        for nm, lvl, mode, obs in manifest:
            n = 0 if mode == "backbone" else K
            fh.write("%s\t%d\t%s\t%s\t%d\n" % (nm, lvl, mode, obs, n))

    print("donor replicon: %s  (len %d)" % (replicon, len(genome)))
    print("native array excised: %d-%d  (%d bp)  ->  backbone %d bp"
          % (ex_lo, ex_hi, ex_hi - ex_lo, len(backbone)))
    print("cassette units used (K): %d   attC per synthetic array: %d" % (K, K))
    print("levels: %s   modes: naive, structure   seed: %d" % (levels, args.seed))
    print("wrote %d genomes to %s/  (+ manifest.tsv, *.truth.tsv)" % (len(manifest), args.out))

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--fna", default="runs_input/vc569B.fna")
    ap.add_argument("--integrons", default="runs/IF_local_vc569B/Results_Integron_Finder_vc569B/vc569B.integrons")
    ap.add_argument("--out", default="synthetic")
    ap.add_argument("--cassettes", type=int, default=80, help="attC per synthetic array")
    ap.add_argument("--levels", default="0,5,10,15,20,25,30", help="target %% divergence levels")
    ap.add_argument("--core", type=int, default=7, help="conserved bp protected at each attC end")
    ap.add_argument("--margin", type=int, default=5000, help="bp excised either side of the native array")
    ap.add_argument("--seed", type=int, default=20260101)
    build(ap.parse_args())
