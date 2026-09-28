"""intervis annotation — name cassette proteins from public databases.

Phase 1: translate each cassette CDS from the genome, DIAMOND-search it against
Swiss-Prot, and take the top hit's protein name. The name is classified into the
same 0-4 function codes the viewer colours by. Each annotator is optional and
fails gracefully, so intervis still renders (as "hypothetical") without them.
"""
from __future__ import annotations
import csv, os, shutil, subprocess, sys, tempfile
from .engine import classify

# minimal standard genetic code
_CODON = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def translate(nt: str) -> str:
    nt = nt.upper().replace("U", "T")
    aa = [_CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3)]
    return "".join(aa).rstrip("*")


def load_fasta(path: str) -> dict:
    """seqid (first whitespace token) -> sequence."""
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def cassette_proteins(parsed: dict, genome_fasta: str) -> dict:
    """element -> protein sequence, translated from the genome by coordinates."""
    seqs = load_fasta(genome_fasta)
    out = {}
    for f in parsed["features"]:
        if f["kind"] != "cassette":
            continue
        rep = seqs.get(f["rep"])
        if rep is None:
            continue
        nt = rep[f["obeg"] - 1:f["oend"]]           # .integrons is 1-based inclusive
        if f["ostrand"] < 0:
            nt = _revcomp(nt)
        aa = translate(nt)
        if aa:
            out[f["element"]] = aa
    return out


# --------------------------------------------------------------------------
# DIAMOND vs Swiss-Prot
# --------------------------------------------------------------------------
def _product(stitle: str) -> str:
    """Pull the readable product description out of a Swiss-Prot subject title."""
    t = stitle.strip().lstrip(">").strip()
    if t.startswith("sp|") or t.startswith("tr|"):        # sp|ACC|NAME desc OS=...
        parts = t.split(None, 1)
        t = parts[1] if len(parts) > 1 else parts[0]
    for cut in (" OS=", " OX=", " GN=", " PE=", " SV="):  # strip UniProt tags
        i = t.find(cut)
        if i != -1:
            t = t[:i]
    return t.strip() or stitle.strip()


def _gene_symbol(stitle: str) -> str:
    """Pull the gene symbol (UniProt GN= field) out of a subject title, or ''."""
    i = stitle.find(" GN=")
    if i == -1:
        return ""
    sym = stitle[i + 4:].split()[0].strip()              # GN=dfrA1 -> dfrA1
    return sym if sym and sym.upper() not in ("NONE", "NA", "-") else ""


def diamond_swissprot(proteins: dict, db: str, threads: int = 4,
                      evalue: float = 1e-5, min_pident: float = 0.0) -> dict:
    """element -> (product, gene, code, pident, evalue) from the best Swiss-Prot hit."""
    if not shutil.which("diamond"):
        raise RuntimeError("diamond not found on PATH — install it "
                           "(e.g. micromamba install -c bioconda diamond)")
    if not (os.path.exists(db) or os.path.exists(db + ".dmnd")):
        raise RuntimeError("Swiss-Prot DIAMOND db not found: %s" % db)
    if not proteins:
        return {}
    tmp = tempfile.mkdtemp(prefix="intervis_")
    qfa, out = os.path.join(tmp, "q.faa"), os.path.join(tmp, "hits.tsv")
    with open(qfa, "w") as fh:
        for el, aa in proteins.items():
            fh.write(">%s\n%s\n" % (el, aa))
    cmd = ["diamond", "blastp", "-q", qfa, "-d", db, "-o", out,
           "--outfmt", "6", "qseqid", "stitle", "pident", "evalue",
           "--max-target-seqs", "1", "-k", "1", "--evalue", str(evalue),
           "--more-sensitive",              # catch divergent cassette homologs
           "-p", str(threads), "--quiet"]
    print("[intervis] annotating %d proteins vs Swiss-Prot..." % len(proteins),
          file=sys.stderr)
    subprocess.run(cmd, check=True)
    res = {}
    with open(out) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 4:
                continue
            el, stitle, pid, ev = row[0], row[1], float(row[2]), float(row[3])
            if el in res or pid < min_pident:
                continue
            product = _product(stitle)
            gene = _gene_symbol(stitle)
            code, _ = classify(product)
            if code == 0:            # a named hit is at least "other"
                code = 1
            res[el] = (product, gene, code, pid, ev)
    shutil.rmtree(tmp, ignore_errors=True)
    return res


# --------------------------------------------------------------------------
# genome annotation overlay (GFF3 / GenBank) — the clinker-style source of
# authentic gene symbols. We read CDS features (seqid, span, strand, gene,
# product) and hand each cassette the annotated CDS it overlaps best.
# --------------------------------------------------------------------------
def _attr(attrs: str, *keys: str) -> str:
    """First matching key= value from a GFF3 attributes column (case-insensitive)."""
    d = {}
    for kv in attrs.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k.strip().lower()] = v.strip()
    for k in keys:
        v = d.get(k.lower())
        if v:
            from urllib.parse import unquote
            return unquote(v)
    return ""


def parse_gff(path: str) -> list:
    """CDS/gene features from a GFF3 -> [(seqid, start, end, strand, gene, product)]."""
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] not in ("CDS", "gene"):
                continue
            try:
                start, end = int(c[3]), int(c[4])
            except ValueError:
                continue
            gene = _attr(c[8], "gene", "Name", "gene_synonym")
            product = _attr(c[8], "product")
            out.append((c[0], min(start, end), max(start, end),
                        c[6], gene, product))
    return out


def parse_genbank(path: str) -> list:
    """CDS features from a GenBank flat file -> [(seqid, start, end, strand, gene, product)].
    A small stdlib parser: enough of the VERSION/FEATURES/CDS grammar to place each
    CDS (feature key at column 6, qualifiers at column 22) and read /gene and /product."""
    import re
    out = []
    seqid = None
    key = None; loc = ""; quals = {}
    qn = None; qv = ""; in_loc = False

    def store_qual():
        nonlocal qn, qv
        if qn:
            quals[qn] = qv.strip().strip('"')
        qn, qv = None, ""

    def emit():
        if key == "CDS" and loc:
            nums = re.findall(r"\d+", loc)
            if nums:
                out.append((seqid, min(int(n) for n in nums), max(int(n) for n in nums),
                            "-" if "complement" in loc else "+",
                            quals.get("gene", ""), quals.get("product", "")))

    with open(path) as fh:
        in_feats = False
        for raw in fh:
            raw = raw.rstrip("\n")
            if raw.startswith("VERSION"):
                p = raw.split(); seqid = p[1] if len(p) > 1 else seqid; continue
            if raw.startswith("ACCESSION") and not seqid:
                p = raw.split(); seqid = p[1] if len(p) > 1 else seqid; continue
            if raw.startswith("FEATURES"):
                in_feats = True; key = None; continue
            if raw.startswith("ORIGIN") or raw.startswith("//"):
                store_qual()
                if key:
                    emit()
                key = None; in_feats = False
                if raw.startswith("//"):
                    seqid = None
                continue
            if not in_feats:
                continue
            if len(raw) >= 6 and raw[:5] == "     " and raw[5] != " ":   # feature key line
                store_qual()
                if key:
                    emit()
                key = raw[5:21].strip(); loc = raw[21:].strip()
                quals = {}; in_loc = True
                continue
            stripped = raw.strip()
            if stripped.startswith("/"):                                 # a qualifier
                store_qual(); in_loc = False
                q = stripped[1:]
                if "=" in q:
                    name, v = q.split("=", 1)
                    if name.strip() in ("gene", "product"):
                        qn, qv = name.strip(), v
                        if len(v) > 1 and v.startswith('"') and v.endswith('"'):
                            store_qual()
            elif qn:                                                     # wrapped qualifier value
                qv += " " + stripped
                if qv.rstrip().endswith('"'):
                    store_qual()
            elif in_loc and key:                                         # wrapped location
                loc += stripped
    return out


def load_annotation(path: str) -> list:
    ext = os.path.splitext(path)[1].lower()
    if ext in (".gff", ".gff3", ".gtf"):
        return parse_gff(path)
    if ext in (".gbff", ".gbk", ".gb", ".gbf", ".genbank"):
        return parse_genbank(path)
    with open(path) as fh:
        head = fh.read(400)
    if head.startswith("LOCUS") or "\nLOCUS" in head:
        return parse_genbank(path)
    return parse_gff(path)                                # default: treat as GFF3


def find_annotation(genome_fasta: str) -> str | None:
    """Look for a matching GFF/GBFF next to the genome FASTA (X.fna -> X.gff...)."""
    stem = os.path.splitext(genome_fasta)[0]
    for ext in (".gff", ".gff3", ".gbff", ".gbk", ".gb"):
        if os.path.exists(stem + ext):
            return stem + ext
    return None


def overlay_map(parsed: dict, anno: list) -> dict:
    """element -> (gene, product) from the annotated CDS each cassette overlaps most."""
    by_rep = {}
    for rec in anno:
        by_rep.setdefault(rec[0], []).append(rec)
    # relaxed replicon lookup: exact id, else ignore a trailing .version
    def recs_for(rep):
        if rep in by_rep:
            return by_rep[rep]
        base = rep.rsplit(".", 1)[0]
        for k, v in by_rep.items():
            if k.rsplit(".", 1)[0] == base:
                return v
        return []
    out = {}
    for f in parsed["features"]:
        if f["kind"] != "cassette":
            continue
        best, bestov = None, 0
        for (_sid, s, e, _st, gene, product) in recs_for(f["rep"]):
            ov = min(f["oend"], e) - max(f["obeg"], s)
            if ov > bestov:
                bestov, best = ov, (gene, product)
        if best and bestov > 0:
            out[f["element"]] = best
    return out


def annotate(parsed: dict, genome_fasta: str, swissprot_db: str | None = None,
             annotation_path: str | None = None, threads: int = 4,
             evalue: float = 1e-5) -> dict:
    """Annotate a parsed genome in place. Gene symbols come from the genome
    annotation (GFF/GBFF) first, then Swiss-Prot's GN= fills the gaps; the
    product drives the function colour. Sets f['gene'], f['product'], f['label'],
    f['code'] on every cassette."""
    sp = {}
    if swissprot_db:
        proteins = cassette_proteins(parsed, genome_fasta)
        sp = diamond_swissprot(proteins, swissprot_db, threads=threads, evalue=evalue)
    anno = {}
    if annotation_path:
        try:
            anno = overlay_map(parsed, load_annotation(annotation_path))
        except Exception as ex:                          # never let a bad file break rendering
            print("[intervis] annotation overlay skipped: %s" % ex, file=sys.stderr)

    named = genes = 0
    for f in parsed["features"]:
        if f["kind"] != "cassette":
            continue
        gene, product = "", None
        if f["element"] in anno:                         # genome annotation wins
            g, p = anno[f["element"]]
            gene = gene or g
            if p and p.lower() not in ("hypothetical protein", "hypothetical"):
                product = product or p
        if f["element"] in sp:                            # Swiss-Prot fills gaps
            sp_product, sp_gene, sp_code, _pid, _ev = sp[f["element"]]
            gene = gene or sp_gene
            product = product or sp_product
        f["gene"] = gene
        f["product"] = product or "hypothetical"
        f["code"], _ = classify(f["product"])
        if gene and f["code"] == 0:                       # a symbol without a known product
            f["code"] = 1
        f["label"] = f["product"]                         # product drives colour + tooltip
        if f["product"] != "hypothetical":
            named += 1
        if gene:
            genes += 1
    n_cass = sum(1 for f in parsed["features"] if f["kind"] == "cassette")
    print("[intervis] %d/%d cassettes with a product, %d with a gene symbol (%d hypothetical)"
          % (named, n_cass, genes, n_cass - named), file=sys.stderr)
    return sp


def write_functions_csv(parsed: dict, path: str) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["element", "gene", "code", "product"])
        for f in parsed["features"]:
            if f["kind"] == "cassette":
                w.writerow([f["element"], f.get("gene", ""),
                            f["code"] if f["code"] is not None else 0,
                            f.get("product") or f["label"] or "hypothetical"])
