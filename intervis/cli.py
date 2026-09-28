"""intervis command-line interface.

  intervis run     GENOME.fna         run IntegronFinder, then render the viewer
  intervis view    FILE.integrons     render the viewer from an existing .integrons
  intervis compare A.integrons B...   render the comparison viewer (2+ arrays)

Phase 0: genome/.integrons -> interactive HTML. Annotation from public
databases and the zoomable attC-in-cassette detail come in later phases.
"""
from __future__ import annotations
import argparse, os, sys
from . import engine
from . import annotate as ann


def _out_path(given: str | None, default_stem: str) -> str:
    if given:
        return given
    return default_stem + ".intervis.html"


def _write(html: str, path: str) -> None:
    with open(path, "w") as fh:
        fh.write(html)
    print("[intervis] wrote", path, "(%d KB)" % (len(html) // 1024), file=sys.stderr)


def cmd_run(a: argparse.Namespace) -> None:
    """One genome -> single-array figure; two (or more) -> clustered comparison.
    Runs IntegronFinder on each genome, then annotates/clusters/renders."""
    if len(a.genome) > 2:
        print("[intervis] note: >2 genomes — fine here; the web tool is limited to 2.",
              file=sys.stderr)
    annos = a.annotation or []
    if annos and len(annos) != len(a.genome):
        sys.exit("--annotation must list one file per genome, in the same order")
    parsed = []
    for i, g in enumerate(a.genome):
        stem = os.path.splitext(os.path.basename(g))[0]
        outdir = os.path.join(a.outdir, stem + "_if") if a.outdir else stem + "_if"
        integrons = engine.run_integron_finder(
            g, outdir, cpu=a.cpu, local_max=not a.no_local_max,
            distance_thresh=a.distance_thresh, evalue_attc=a.evalue_attc)
        p = engine.parse_integrons(integrons, bin_id=(a.name if len(a.genome) == 1 else stem))
        _report(p)
        engine.attach_sequence(p, g)              # carry the region sequence for base-level zoom
        if annos:                                 # explicit list ("NONE" = skip this genome)
            anno = annos[i] if annos[i] not in ("", "NONE", "-") else None
        else:
            anno = ann.find_annotation(g)         # else a matching GFF/GBFF beside the FASTA
        if a.annotate or anno:
            ann.annotate(p, g, swissprot_db=(a.swissprot if a.annotate else None),
                         annotation_path=anno, threads=a.cpu)
        parsed.append(p)

    out = _out_path(a.output, parsed[0]["bin_id"] if len(parsed) == 1 else "comparison")
    if len(parsed) == 1:
        if a.annotate:
            ann.write_functions_csv(parsed[0], os.path.splitext(out)[0] + ".functions.csv")
        _write(engine.build_cartographer(parsed[0]), out)
    else:
        from . import cluster
        fams, idn, famfn, links = cluster.cluster_families(parsed, a.genome, threads=a.cpu)
        _write(engine.build_synteny(parsed, fams, idn, famfn, links), out)


def cmd_view(a: argparse.Namespace) -> None:
    parsed = engine.parse_integrons(a.integrons, bin_id=a.name)
    _report(parsed)
    out = _out_path(a.output, parsed["bin_id"])
    anno = a.annotation
    if a.genome:
        engine.attach_sequence(parsed, a.genome)  # base-level zoom needs the genome
        if not anno:
            anno = ann.find_annotation(a.genome)  # a matching GFF/GBFF next to the FASTA
    if a.annotate or anno:
        if a.annotate and not a.genome:
            sys.exit("--annotate needs --genome (the FASTA the .integrons came from)")
        if anno and not a.genome:
            print("[intervis] using genome annotation for gene names", file=sys.stderr)
        ann.annotate(parsed, a.genome, swissprot_db=(a.swissprot if a.annotate else None),
                     annotation_path=anno, threads=a.cpu)
        ann.write_functions_csv(parsed, os.path.splitext(out)[0] + ".functions.csv")
    _write(engine.build_cartographer(parsed), out)


def cmd_compare(a: argparse.Namespace) -> None:
    if len(a.integrons) < 2:
        sys.exit("compare needs at least two .integrons files")
    if len(a.integrons) > 2:
        print("[intervis] note: >2 arrays — fine on the command line; the web "
              "tool is limited to 2.", file=sys.stderr)
    parsed = [engine.parse_integrons(p) for p in a.integrons]
    for p in parsed:
        _report(p)
    families = identity = fam_fn = links = None
    if a.genomes:
        if len(a.genomes) != len(a.integrons):
            sys.exit("--genomes must list one FASTA per .integrons, in the same order")
        annos = a.annotations or []
        if annos and len(annos) != len(a.genomes):
            sys.exit("--annotations must list one file per genome, in the same order")
        for i, (p, g) in enumerate(zip(parsed, a.genomes)):
            engine.attach_sequence(p, g)      # carry each region's sequence for base-level zoom
            if annos:
                anno = annos[i] if annos[i] not in ("", "NONE", "-") else None
            else:
                anno = ann.find_annotation(g)
            if a.swissprot or anno:           # name the families too
                ann.annotate(p, g, swissprot_db=a.swissprot, annotation_path=anno, threads=a.cpu)
        from . import cluster
        families, identity, fam_fn, links = cluster.cluster_families(
            parsed, a.genomes, threads=a.cpu, min_id=a.min_id, min_cov=a.min_cov)
    else:
        print("[intervis] note: no --genomes given, so ribbons need --genomes to link "
              "cassettes by homology. Pass --genomes A.fna B.fna for real ribbons.",
              file=sys.stderr)
    _write(engine.build_synteny(parsed, families, identity, fam_fn, links),
           _out_path(a.output, "comparison"))


def _report(parsed: dict) -> None:
    feats = parsed["features"]
    from collections import Counter
    kinds = Counter(f["kind"] for f in feats)
    print("[intervis] %s: %d cassettes, %d attC, span %.1f kb"
          % (parsed["bin_id"], kinds.get("cassette", 0), kinds.get("attC", 0),
             parsed["length"] / 1000.0), file=sys.stderr)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="intervis",
                                description="Interactive integron array visualiser.")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run IntegronFinder on 1-2 genomes, then render")
    r.add_argument("genome", nargs="+", help="one genome FASTA -> single view; two -> comparison")
    r.add_argument("-o", "--output", help="output HTML path")
    r.add_argument("--outdir", help="IntegronFinder output dir (default <genome>_if)")
    r.add_argument("--name", help="label shown in the viewer (default: file stem)")
    r.add_argument("--cpu", type=int, default=4)
    r.add_argument("--no-local-max", action="store_true",
                   help="disable --local-max (not recommended; lowers sensitivity)")
    r.add_argument("--distance-thresh", type=int, default=None)
    r.add_argument("--evalue-attc", type=float, default=None)
    r.add_argument("--annotate", action="store_true",
                   help="name cassettes from Swiss-Prot (needs --swissprot, diamond on PATH)")
    r.add_argument("--swissprot", help="path to the Swiss-Prot DIAMOND db (.dmnd)")
    r.add_argument("--annotation", nargs="+",
                   help="genome annotation (GFF3/GenBank) per genome, same order — real "
                        "gene symbols by coordinate overlap (else a matching file is auto-found)")
    r.set_defaults(func=cmd_run)

    v = sub.add_parser("view", help="render from an existing .integrons file")
    v.add_argument("integrons", help="path to a .integrons file")
    v.add_argument("-o", "--output", help="output HTML path")
    v.add_argument("--name", help="label shown in the viewer (default: file stem)")
    v.add_argument("--annotate", action="store_true",
                   help="name cassettes from Swiss-Prot (needs --genome, --swissprot, diamond)")
    v.add_argument("--genome", help="genome FASTA the .integrons came from (for --annotate)")
    v.add_argument("--swissprot", help="path to the Swiss-Prot DIAMOND db (.dmnd)")
    v.add_argument("--annotation", help="genome annotation (GFF3/GenBank) for gene names")
    v.add_argument("--cpu", type=int, default=4)
    v.set_defaults(func=cmd_view)

    c = sub.add_parser("compare", help="render a comparison of 2+ arrays")
    c.add_argument("integrons", nargs="+", help="two or more .integrons files")
    c.add_argument("-o", "--output", help="output HTML path")
    c.add_argument("--genomes", nargs="+",
                   help="one genome FASTA per .integrons, same order — enables real "
                        "cross-genome family clustering (ribbons)")
    c.add_argument("--swissprot", help="Swiss-Prot DIAMOND db, to also name the families")
    c.add_argument("--annotations", nargs="+",
                   help="genome annotation (GFF3/GenBank) per genome, same order — gene names")
    c.add_argument("--cpu", type=int, default=4)
    c.add_argument("--min-id", type=float, default=30.0,
                   help="min %% identity to link cassettes into a family (default 30)")
    c.add_argument("--min-cov", type=float, default=0.5,
                   help="min alignment coverage of both proteins (default 0.5)")
    c.set_defaults(func=cmd_compare)
    return p


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
