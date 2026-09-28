"""intervis cross-genome cassette clustering.

Groups the cassette proteins of two (or more) genomes into families by an
all-vs-all DIAMOND search + single-linkage (union-find) clustering. Families
shared between genomes are what the comparison viewer draws ribbons for, shaded
by mean cross-genome % identity. Reuses the DIAMOND already installed for
annotation — no extra tool. Falls back cleanly (raises) if DIAMOND is absent.
"""
from __future__ import annotations
import csv, os, shutil, subprocess, sys, tempfile
from collections import defaultdict
from .annotate import cassette_proteins


class _UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:      # path compression
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def _gather(parsed_list, genome_paths):
    """(genome_index, element) proteins, keyed 'gi::element' for DIAMOND."""
    prots, meta = {}, {}
    for gi, (parsed, gpath) in enumerate(zip(parsed_list, genome_paths)):
        for el, aa in cassette_proteins(parsed, gpath).items():
            key = "%d::%s" % (gi, el)
            prots[key] = aa
            meta[key] = (gi, el)
    return prots, meta


def cluster_families(parsed_list, genome_paths, threads: int = 4,
                     min_id: float = 30.0, min_cov: float = 0.5) -> tuple:
    """Return (families, identity, fam_fn, links):
        families : {(gi, element): famid}
        identity : {famid: 0..1}     mean cross-genome identity (shared fams only)
        fam_fn   : {famid: [code, label]}
        links    : [(gi, el_a, gj, el_b, pident, cov)]  reciprocal-best-hit cassette
                   pairs between genomes — each with its EXACT pairwise % identity and
                   coverage (this is what the ribbons are drawn and coloured from, so a
                   genome vs itself gives 100% for every cassette's identical twin).
    """
    prots, meta = _gather(parsed_list, genome_paths)
    if not prots:
        return {}, {}, {}
    if not shutil.which("diamond"):
        raise RuntimeError("diamond not found on PATH (needed for clustering)")

    tmp = tempfile.mkdtemp(prefix="intervis_clu_")
    faa, db, out = (os.path.join(tmp, x) for x in ("all.faa", "all", "hits.tsv"))
    with open(faa, "w") as fh:
        for k, aa in prots.items():
            fh.write(">%s\n%s\n" % (k, aa))
    print("[intervis] clustering %d cassette proteins across %d genomes..."
          % (len(prots), len(parsed_list)), file=sys.stderr)
    subprocess.run(["diamond", "makedb", "--in", faa, "-d", db, "--quiet"], check=True)
    subprocess.run(["diamond", "blastp", "-q", faa, "-d", db, "-o", out, "--outfmt", "6",
                    "qseqid", "sseqid", "pident", "length", "qlen", "slen", "bitscore",
                    "--more-sensitive", "-k", "50", "--evalue", "1e-5",
                    "-p", str(threads), "--quiet"], check=True)

    uf = _UnionFind()
    for k in prots:
        uf.find(k)                                   # every protein is at least a singleton
    edges = []                                       # (q, s, pident) meeting thresholds
    best = {}                                        # (query, target_genome) -> best cross-genome hit
    with open(out) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 7:
                continue
            q, s = row[0], row[1]
            pid, ln, ql, sl = float(row[2]), int(row[3]), int(row[4]), int(row[5])
            bits = float(row[6])
            if q == s:
                continue
            gi_q, gi_s = meta[q][0], meta[s][0]
            cov = ln / min(ql, sl)                    # alignment coverage of the shorter protein
            if pid >= min_id and ln / ql >= min_cov and ln / sl >= min_cov:
                uf.union(q, s)
                if gi_q != gi_s:
                    edges.append((q, s, pid))
            # track each query's best hit into each OTHER genome (by bitscore) for RBH
            if gi_q != gi_s and pid >= min_id and cov >= min_cov:
                key = (q, gi_s)
                if key not in best or bits > best[key][3]:
                    best[key] = (s, pid, cov, bits)

    # reciprocal best hits -> the orthologous cassette pairs the ribbons connect,
    # each carrying its OWN pairwise % identity (not a family average)
    links, seen = [], set()
    for (q, gj), (s, pid, cov, _bits) in best.items():
        gi = meta[q][0]
        back = best.get((s, gi))
        if not back or back[0] != q:                 # must be mutual best hits
            continue
        pk = tuple(sorted((q, s)))
        if pk in seen:
            continue
        seen.add(pk)
        (giq, elq), (gis, els) = meta[q], meta[s]
        links.append((giq, elq, gis, els, round(pid, 1), round(cov, 3)))

    # contiguous family ids
    root_fam, families = {}, {}
    for key, (gi, el) in meta.items():
        fid = root_fam.setdefault(uf.find(key), len(root_fam))
        families[(gi, el)] = fid

    # mean cross-genome identity per (shared) family -> ribbon shade
    fam_edges = defaultdict(list)
    for q, s, pid in edges:
        if meta[q][0] != meta[s][0] and families[meta[q]] == families[meta[s]]:
            fam_edges[families[meta[q]]].append(pid)
    identity = {f: round(sum(v) / len(v) / 100.0, 3) for f, v in fam_edges.items() if v}

    # family function = the best-annotated member (prefer a real function over hypothetical)
    ann = {}
    for gi, parsed in enumerate(parsed_list):
        for f in parsed["features"]:
            if f["kind"] == "cassette":
                ann[(gi, f["element"])] = [f["code"] or 0, f["label"] or "hypothetical"]
    fam_fn = {}
    for gikey, fid in families.items():
        cl = ann.get(gikey, [0, "hypothetical"])
        if fid not in fam_fn or (fam_fn[fid][0] == 0 and cl[0] != 0):
            fam_fn[fid] = cl

    shutil.rmtree(tmp, ignore_errors=True)
    n_shared = sum(1 for f in set(families.values())
                   if len({gi for (gi, _), ff in families.items() if ff == f}) >= 2)
    print("[intervis] %d families, %d shared across genomes; %d reciprocal-best-hit links"
          % (len(set(families.values())), n_shared, len(links)), file=sys.stderr)
    return families, identity, fam_fn, links
