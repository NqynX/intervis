"""intervis engine — IntegronFinder output -> interactive viewer.

Phase 0: run IntegronFinder on a genome (or take an existing .integrons),
parse it into oriented features, and render the interactive HTML viewer.
Stdlib only. The parse is a faithful port of the validated pipeline logic
(parse.R / augment_comparison.py): integrase-left orientation, "pick the
array with the most attC", kinds integrase / cassette / attC.
"""
from __future__ import annotations
import csv, glob, json, os, re, subprocess, sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATES = os.path.join(HERE, "templates")

# function code -> palette slot (matches the viewer legend)
#   0 hypothetical  1 other/Pfam  2 toxin-antitoxin  3 anti-phage defence  4 AMR
TA_KEYS = ("toxin", "antitoxin", "addiction", "pare", "pard", "phd", "yefm",
           "rele", "relb", "vapc", "vapb", "maze", "mazf", "hica", "hicb",
           "ccda", "ccdb", "higa", "higb")
AMR_KEYS = ("lactamas", "aminoglycoside", "chloramphenicol", "dihydrofolate",
            "dihydropteroate", "quinolone", "qnr", "tetracycline", "sulfona",
            "rifampin", "macrolide", "dfra", "aac(", "aph(", "resistance")
DEF_KEYS = ("defense", "defence", "anti-phage", "antiphage", "restriction",
            "abortive", "crispr", "cas", "retron", "gabija", "hachiman")


def classify(text: str) -> tuple[int, str]:
    """Map an annotation string to (code, label). Phase-0 keyword pass."""
    t = (text or "").strip()
    low = t.lower()
    if not low or low in ("protein", "hypothetical protein", "na", "-"):
        return 0, "hypothetical"
    if any(k in low for k in DEF_KEYS):
        return 3, t
    if any(k in low for k in TA_KEYS):
        return 2, t
    if any(k in low for k in AMR_KEYS):
        return 4, t
    return 1, t


# --------------------------------------------------------------------------
# parse
# --------------------------------------------------------------------------
NEED = ("ID_integron", "ID_replicon", "element", "pos_beg", "pos_end",
        "strand", "type_elt", "annotation", "model")


def parse_integrons(path: str, integron: str | None = None,
                    bin_id: str | None = None) -> dict:
    """Parse one integron from a `.integrons` file into oriented features.

    Returns {bin_id, length, features:[...]} where each feature is
    {element, kind, start, end, strand, code, label} with kind in
    {integrase, cassette, attC}, oriented integrase-left, sorted by start.
    """
    if bin_id is None:
        bin_id = os.path.splitext(os.path.basename(path))[0]
    rows, hdr = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = f
                continue
            rows.append(dict(zip(hdr, f)))
    if hdr is None or any(c not in hdr for c in NEED):
        miss = [c for c in NEED if hdr is None or c not in hdr]
        raise ValueError("Not a valid .integrons file (missing columns: %s)"
                         % ", ".join(miss))
    for r in rows:
        r["pos_beg"] = int(float(r["pos_beg"]))
        r["pos_end"] = int(float(r["pos_end"]))
        r["key"] = r["ID_replicon"] + "|" + r["ID_integron"]

    # pick the integron with the most attC sites (fall back to any)
    if integron is None:
        attc_keys = [r["key"] for r in rows if r["type_elt"].lower() == "attc"]
        tally = Counter(attc_keys) or Counter(r["key"] for r in rows)
        integron = tally.most_common(1)[0][0]
    sel = [r for r in rows if r["key"] == integron]
    if not sel:
        raise ValueError("Integron %r not found in %s" % (integron, path))

    a0 = min(r["pos_beg"] for r in sel)
    a1 = max(r["pos_end"] for r in sel)
    isinti = lambda r: "inti" in (r["annotation"] + " " + r["model"]).lower()
    inti = [r for r in sel if isinti(r)]
    flip = bool(inti) and (sum(r["pos_beg"] for r in inti) / len(inti) > (a0 + a1) / 2)

    feats = []
    for r in sel:
        s, e = r["pos_beg"], r["pos_end"]
        d = -1 if r["strand"] in ("-1", "-") else 1
        if flip:
            s, e = a1 - e + a0, a1 - s + a0
            d = -d
        if r["type_elt"].lower() == "attc":
            kind, code, label = "attC", None, None
        elif isinti(r):
            kind, code, label = "integrase", None, None
        else:
            kind = "cassette"
            code, label = classify(r.get("annotation", ""))
        feats.append({"element": r["element"], "kind": kind,
                      "start": min(s, e), "end": max(s, e),
                      "strand": "-" if d < 0 else "+",
                      "code": code, "label": label,
                      # original (un-flipped) coordinates, for protein extraction
                      "rep": r["ID_replicon"], "obeg": r["pos_beg"],
                      "oend": r["pos_end"],
                      "ostrand": -1 if r["strand"] in ("-1", "-") else 1})
    feats.sort(key=lambda f: f["start"])
    return {"bin_id": bin_id, "length": a1 - a0, "features": feats,
            "a0": a0, "a1": a1, "flip": flip, "rep": sel[0]["ID_replicon"]}


_NT_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(s: str) -> str:
    return s.translate(_NT_COMP)[::-1]


def _load_fasta(path: str) -> dict:
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def attach_sequence(parsed: dict, genome_fasta: str) -> None:
    """Attach the integron region's nucleotide sequence, oriented to match the
    displayed (integrase-left) coordinates, so the viewer can show bases on zoom.
    Sets parsed['seq'] (string) and parsed['seq0'] (its genomic start = a0)."""
    rep = _load_fasta(genome_fasta).get(parsed["rep"])
    if rep is None:
        return
    a0, a1 = parsed["a0"], parsed["a1"]
    region = rep[a0 - 1:a1].upper()             # .integrons is 1-based inclusive
    if parsed["flip"]:
        region = _revcomp(region)               # match the flipped display frame
    parsed["seq"] = region
    parsed["seq0"] = a0


def find_integrons(results_dir: str) -> str:
    """Locate the .integrons file inside an IntegronFinder --outdir."""
    hits = glob.glob(os.path.join(results_dir, "**", "*.integrons"),
                     recursive=True)
    if not hits:
        raise FileNotFoundError("no .integrons under %s" % results_dir)
    return hits[0]


# --------------------------------------------------------------------------
# run IntegronFinder
# --------------------------------------------------------------------------
def run_integron_finder(genome: str, outdir: str, cpu: int = 4,
                        local_max: bool = True, keep_palindromes: bool = True,
                        distance_thresh: int | None = None,
                        evalue_attc: float | None = None,
                        extra: list[str] | None = None) -> str:
    """Run IntegronFinder on a genome FASTA; return the .integrons path.

    Uses the same flags as the validated pipeline. Requires `integron_finder`
    on PATH (activate the integronfinder env first).
    """
    os.makedirs(outdir, exist_ok=True)
    cmd = ["integron_finder", genome, "--outdir", outdir, "--cpu", str(cpu)]
    if local_max:
        cmd.append("--local-max")
    if keep_palindromes:
        cmd.append("--keep-palindromes")
    if distance_thresh is not None:
        cmd += ["--distance-thresh", str(distance_thresh)]
    if evalue_attc is not None:
        cmd += ["--evalue-attc", str(evalue_attc)]
    if extra:
        cmd += extra
    print("[intervis] running:", " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)
    return find_integrons(outdir)


# --------------------------------------------------------------------------
# render
# --------------------------------------------------------------------------
KIND_CODE = {"integrase": 0, "cassette": 1, "attC": 2}


def _inject(html: str, marker: str, value: str) -> str:
    i = html.find(marker)
    if i < 0:
        raise ValueError("placeholder %s not found in template" % marker)
    return html[:i] + value + html[i + len(marker):]


def _read_template(name: str) -> str:
    with open(os.path.join(TEMPLATES, name)) as fh:
        return fh.read()


def build_cartographer(parsed: dict) -> str:
    """Single array -> Integron Cartographer HTML."""
    feats = parsed["features"]
    FEAT = [[f["start"], f["end"], -1 if f["strand"] == "-" else 1,
             KIND_CODE[f["kind"]]] for f in feats]
    FN = [[f["code"] if f["code"] is not None else 0,
           f.get("product") or f["label"] or "hypothetical",
           f.get("gene", "")]
          for f in feats if f["kind"] == "cassette"]
    html = _read_template("cartographer.html")
    html = _inject(html, "__DATA__", json.dumps(FEAT))
    html = _inject(html, "__FN__", json.dumps(FN))
    html = _inject(html, "__NAME__", json.dumps(parsed["bin_id"]))
    html = _inject(html, "__SEQ__", json.dumps(parsed.get("seq", "")))
    html = _inject(html, "__SEQ0__", json.dumps(parsed.get("seq0", 0)))
    return html


def build_synteny(parsed_list: list[dict], families: dict | None = None,
                  identity: dict | None = None, fam_fn: dict | None = None,
                  links: list | None = None) -> str:
    """Two (or more) arrays -> Integron Synteny HTML.

    families : {(genome_index, element): famid} from cross-genome clustering
               (cluster.py). When None, falls back to keying by identical
               element id (only meaningful for re-runs of the same genome).
    identity : {famid: 0..1} mean cross-genome % identity -> ribbon shade.
    fam_fn   : {famid: [code, label]} function/name per family.
    links    : [(gi, el_a, gj, el_b, pident, cov)] reciprocal-best-hit cassette
               pairs -> the ribbons, each with its exact pairwise % identity.
    """
    fam_ids: dict[str, int] = {}
    genomes = []
    fam_dict: dict[str, list] = {}
    for gi, parsed in enumerate(parsed_list):
        rows = []
        for f in parsed["features"]:
            kc = KIND_CODE[f["kind"]]
            fid = -1
            if f["kind"] == "cassette":
                if families is not None:
                    fid = families.get((gi, f["element"]), -1)
                else:
                    fid = fam_ids.setdefault(f["element"], len(fam_ids))
                if fid != -1:
                    if fam_fn and fid in fam_fn:
                        fam_dict.setdefault(str(fid), fam_fn[fid])
                    else:
                        fam_dict.setdefault(str(fid),
                                            [f["code"] or 0, f["label"] or "hypothetical"])
            code = (f["code"] if f.get("code") is not None else 0)
            label = f.get("product") or f.get("label") or ""
            gene = f.get("gene", "") or ("intI" if f["kind"] == "integrase" else "")
            rows.append([f["start"], f["end"],
                         -1 if f["strand"] == "-" else 1, kc, fid,
                         code, label, gene])       # + per-cassette function/name (cartographer style)
        rows.sort(key=lambda z: z[0])
        genomes.append({"name": parsed["bin_id"], "acc": parsed["bin_id"], "f": rows,
                        "seq": parsed.get("seq", ""), "seq0": parsed.get("seq0", 0)})
    pid = {str(k): v for k, v in (identity or {}).items()}
    # map each cassette element to its display span (integrase-left coords) per genome,
    # then translate the reciprocal-best-hit links into drawable ribbon endpoints
    elpos = []
    for parsed in parsed_list:
        m = {}
        for f in parsed["features"]:
            if f["kind"] == "cassette":
                m[f["element"]] = (f["start"], f["end"])
        elpos.append(m)
    link_rows = []
    for (gi, ela, gj, elb, pident, cov) in (links or []):
        a = elpos[gi].get(ela) if 0 <= gi < len(elpos) else None
        b = elpos[gj].get(elb) if 0 <= gj < len(elpos) else None
        if not a or not b:
            continue
        fam = (families or {}).get((gi, ela), -1)
        link_rows.append([gi, a[0], a[1], gj, b[0], b[1], pident, cov, fam])
    comp = {"genomes": genomes, "fam": fam_dict, "pid": pid, "links": link_rows}
    html = _read_template("synteny.html")
    return _inject(html, "__COMP__", json.dumps(comp))
