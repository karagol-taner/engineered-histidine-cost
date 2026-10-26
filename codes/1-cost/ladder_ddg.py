#!/usr/bin/env python3
"""Cost of each substitution against a single-nucleotide-reachable baseline.

Reads the cDNA-proteolysis, FireProtDB and S669 stability sets, the AlphaMissense
cohorts and the UniProt humsavar variant file, maps every stability column to a
common sign, and writes the per-substitution Cliff delta for each readout.

    python3 ladder_ddg.py --help
"""
import argparse, bisect, csv, glob, io, importlib.util, os, re, statistics as st, subprocess, sys

MUT = re.compile(r'^([A-Z])(-?\d+)([A-Z])$')
AA = set("ACDEFGHIKLMNPQRSTVWY")


def tpath(targets, filename):
    """Where --targets would keep one of the method's own modules, else ""."""
    return os.path.join(targets, filename) if targets else ""


def load(path, name):
    """The substitution set and the single-nucleotide reachability test.

    Reads the bundled codes/ladder_sets.py, which carries the substitution set
    and computes single-nucleotide reachability from the standard genetic
    code. --targets overrides it with the method's own modules,
    which are not distributed here; the two give identical answers on every
    name used, checked over all 400 ordered amino-acid pairs for reachability.
    """
    if not path or not os.path.exists(path):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            os.pardir, "ladder_sets.py")
        name = "ladder_sets"
        if not os.path.exists(path):
            sys.exit("not found: %s" % path)
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ------------------------------------------------------------------- readers
def read_fireprot(d):
    """(aa1, aa2, raw_ddG, asa, pH) with duplicates removed."""
    seen, out = set(), []
    for p in sorted(glob.glob(os.path.join(d, "*.csv"))):
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh):
                a1, a2, v = r.get("wild_type"), r.get("mutation"), r.get("ddG")
                if not a1 or not a2 or a1 not in AA or a2 not in AA or a1 == a2:
                    continue
                try:
                    v = float(v)
                except (TypeError, ValueError):
                    continue
                key = (r.get("uniprot_id"), r.get("position"), a1, a2, round(v, 3))
                if key in seen:
                    continue
                seen.add(key)
                try:
                    asa = float(r.get("asa") or "nan")
                except ValueError:
                    asa = float("nan")
                try:
                    ph = float(r.get("pH") or "nan")
                except ValueError:
                    ph = float("nan")
                out.append((a1, a2, v, asa, ph))
    return out


def read_megascale(d):
    out = []
    for p in sorted(glob.glob(os.path.join(d, "*.csv"))):
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh):
                m = MUT.match((r.get("mut_type") or "").strip())
                if not m:
                    continue
                a1, _, a2 = m.groups()
                if a1 not in AA or a2 not in AA or a1 == a2:
                    continue
                v = r.get("ddG_ML")
                try:
                    v = float(v)
                except (TypeError, ValueError):
                    continue
                out.append((a1, a2, v, float("nan"), float("nan")))
    return out


def read_s669(d):
    out = []
    for p in sorted(glob.glob(os.path.join(d, "*.csv"))):
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh):
                m = MUT.match((r.get("MUT") or "").strip())
                if not m:
                    continue
                a1, _, a2 = m.groups()
                if a1 not in AA or a2 not in AA or a1 == a2:
                    continue
                try:
                    v = float(r.get("DDG"))
                except (TypeError, ValueError):
                    continue
                out.append((a1, a2, v, float("nan"), float("nan")))
    return out


# -------------------------------------------------------------- sign mapping
MIN_PRO = 30          # below this the proline check is noise, not a check


def class_means(cost, minn=3):
    per = {}
    for a1, a2, v, _, _ in cost:
        per.setdefault((a1, a2), []).append(v)
    return {k: st.mean(v) for k, v in per.items() if len(v) >= minn}


def to_cost(rows, tag, reference=None):
    """Map each source's own sign to one convention: higher = more destabilising."""
    vals = [r[2] for r in rows]
    frac_neg = sum(1 for v in vals if v < 0) / len(vals)
    flip = frac_neg > 0.5
    cost = [(a1, a2, (-v if flip else v), asa, ph) for a1, a2, v, asa, ph in rows]
    print("%s: n=%d, %.1f%% of values negative, so %s is destabilising; "
          "mapped to cost with higher = less stable"
          % (tag, len(rows), 100 * frac_neg, "negative" if flip else "positive"))

    med = st.median([c[2] for c in cost])
    pro = [c[2] for c in cost if c[1] == "P"]
    if len(pro) >= MIN_PRO:
        mpro = st.median(pro)
        print("   sign check, proline: median cost %+.3f overall against %+.3f for "
              "X>P (n=%d)" % (med, mpro, len(pro)))
        if mpro <= med:
            sys.exit("FATAL %s: substitutions to proline are not costlier than the "
                     "median, so the sign mapping is inverted or the column is not "
                     "a stability change. Refusing to rank edges on it." % tag)
        return cost

    # Too few prolines for that test.
    if not reference:
        print("   only %d X>P substitutions and no reference set, so the sign is "
              "unchecked; this source is reported and not relied on" % len(pro))
        return cost
    mine, theirs = class_means(cost), class_means(reference)
    shared = sorted(set(mine) & set(theirs))
    if len(shared) < 20:
        print("   only %d substitution classes shared with the reference, sign "
              "unchecked" % len(shared))
        return cost
    rho = spearman([mine[k] for k in shared], [theirs[k] for k in shared])
    print("   sign check, against the reference set: rho = %+.3f over %d shared "
          "substitution classes (%d X>P, too few for the proline test)"
          % (rho, len(shared), len(pro)))
    if rho < 0.3:
        sys.exit("FATAL %s: per-class costs do not agree with the reference set "
                 "(rho %+.3f), so the sign mapping or the column is wrong." % (tag, rho))
    return cost


# -------------------------------------------------------------------- scoring
def cliff_factory(base):
    base = sorted(base)
    nb = len(base)

    def cliff(xs):
        gt = lt = 0
        for x in xs:
            lt += bisect.bisect_left(base, x)
            gt += nb - bisect.bisect_right(base, x)
        return (gt - lt) / (len(xs) * nb)
    return cliff, nb, base


def spearman(x, y):
    def rank(v):
        s = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(s):
            j = i
            while j + 1 < len(s) and v[s[j + 1]] == v[s[i]]:
                j += 1
            for k in range(i, j + 1):
                r[s[k]] = (i + j) / 2.0 + 1
            i = j + 1
        return r
    rx, ry = rank(x), rank(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else float("nan")


# ------------------------------------------------- the two published rankings
def run_and_parse(script, datafile, outpath, edge_i, val_i, label):
    """Run one of the method's own scripts, or read its deposited output.

    With --targets set the script is run and its raw output is written to
    outpath. Without it the deposited outpath is read instead, so the
    comparison reproduces from this repository alone.
    """
    out = None
    if script and os.path.exists(script):
        try:
            out = subprocess.run([sys.executable, script, datafile], check=True,
                                 capture_output=True, text=True,
                                 timeout=1800).stdout
            with open(outpath, "w") as fh:
                fh.write(out)
        except Exception as e:
            print("   could not run %s: %s" % (os.path.basename(script), e))
            out = None
    if out is None:
        if not os.path.exists(outpath):
            print("   no deposited output at %s, and --targets is not set" % outpath)
            return {}
        print("   read the deposited output of %s; set --targets to recompute it"
              % os.path.basename(outpath))
        out = io.open(outpath, encoding="utf-8").read()
    vals, inblock = {}, False
    for line in out.splitlines():
        if line.startswith("--- 12 core edges"):
            inblock = True; continue
        if inblock and line.startswith("---"):
            break
        f = line.split()
        if inblock and len(f) > max(edge_i, val_i) and ">" in f[edge_i]:
            try:
                vals[tuple(f[edge_i].split(">"))] = float(f[val_i])
            except ValueError:
                pass
    if len(vals) != 12:
        sys.exit("FATAL: parsed %d of 12 core edges from %s. The output format "
                 "changed, and a correlation computed on a partial table would be "
                 "wrong rather than merely incomplete." % (len(vals), label))
    print("   %s: 12 core edges parsed, raw output kept at %s"
          % (label, os.path.basename(outpath)))
    return vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fireprot", default="../../data/fireprot")
    ap.add_argument("--megascale", default="../../data/megascale")
    ap.add_argument("--s669", default="../../data/s669")
    ap.add_argument("--targets", default="",
                    help="directory holding the method's own "
                         "modules; omit to use the bundled "
                         "codes/ladder_sets.py")
    ap.add_argument("--humsavar", default="../../data/humsavar.txt")
    ap.add_argument("--alphamissense", default="../../data/alphamissense_soluble.csv.gz")
    ap.add_argument("--outdir", default="../../data")
    ap.add_argument("--buried-asa", type=float, default=25.0)
    a = ap.parse_args()

    lec = load(tpath(a.targets, "ladder_edge_cost.py"), "lec")
    ca = load(tpath(a.targets, "codon_audit.py"), "ca")
    CORE, ARG, EXT = lec.CORE, lec.ARG, lec.EXT
    BASICDIR, CHARGE = lec.BASICDIR, lec.CHARGE
    snv_ok = lambda p: bool(ca.routes(p[0], p[1]))

    print("ladder_ddg.py\n")

    sets, ref = [], None
    for tag, reader, d in (("megascale", read_megascale, a.megascale),
                           ("fireprot", read_fireprot, a.fireprot),
                           ("s669", read_s669, a.s669)):
        rows = reader(d)
        if not rows:
            print("%s: no rows read from %s, skipped" % (tag, d)); continue
        c = to_cost(rows, tag, reference=ref)
        if ref is None or len(c) > len(ref):
            ref = c
        sets.append((tag, c))
    if not sets:
        sys.exit("no measured data read at all")
    print()

    results = {}
    for tag, cost in sets:
        base = [c[2] for c in cost if snv_ok((c[0], c[1]))]
        if not base:
            print("%s: no single-nucleotide-reachable substitutions, skipped" % tag); continue
        cliff, nb, _ = cliff_factory(base)
        per = {}
        for a1, a2, v, asa, ph in cost:
            per.setdefault((a1, a2), []).append((v, asa))
        print("=" * 100)
        print("%s: baseline is every single-nucleotide-reachable substitution, "
              "n=%d, mean %+.3f kcal/mol, median %+.3f"
              % (tag, nb, st.mean(base), st.median(base)))
        print("positive Cliff delta = less destabilising than an average "
              "single-base substitution")
        hdr = ("%-7s %-7s %-7s %-9s %-9s %-8s %-7s %s"
               % ("edge", "dir", "n", "mean", "median", "cliff", "dQ", "buried / exposed mean"))
        for pairs, name in ((CORE, "12 core edges, best tolerated first"),
                            (ARG, "5 arginine edges"),
                            (EXT, "proposed extensions")):
            print("\n--- %s ---" % name); print(hdr)
            rows_out = []
            for p in pairs:
                v = per.get(p)
                if not v or len(v) < 3:
                    continue
                vals = [x for x, _ in v]
                d = cliff(vals)
                bu = [x for x, s in v if s == s and s < a.buried_asa]
                ex = [x for x, s in v if s == s and s >= a.buried_asa]
                strat = ("%+.2f (n=%d) / %+.2f (n=%d)"
                         % (st.mean(bu), len(bu), st.mean(ex), len(ex))) if bu and ex else "-"
                rows_out.append((-d, "%-7s %-7s %-7d %+-9.3f %+-9.3f %+-8.3f %+-7.1f %s"
                                 % (p[0] + ">" + p[1],
                                    "basic" if p in BASICDIR else
                                    "acidic" if (p[1], p[0]) in BASICDIR else "-",
                                    len(vals), st.mean(vals), st.median(vals), d,
                                    CHARGE[p[1]] - CHARGE[p[0]], strat)))
                if pairs is CORE:
                    results.setdefault(tag, {})[p] = d
            for _, l in sorted(rows_out):
                print(l)
        print()

    # ------------------------------------------------- against the two rankings
    print("=" * 100)
    print("agreement with the two published rankings, both recomputed from source")
    am = run_and_parse(tpath(a.targets, "ladder_edge_cost.py"), a.alphamissense,
                       os.path.join(a.outdir, "source_runs", "ladder_edge_cost_soluble.out"),
                       0, 6, "AlphaMissense soluble, Cliff delta")
    hs = run_and_parse(tpath(a.targets, "ladder_humsavar.py"), a.humsavar,
                       os.path.join(a.outdir, "source_runs", "ladder_humsavar.out"),
                       0, 5, "humsavar, pathogenic per cent")
    print()
    for tag, edges in results.items():
        common = [p for p in CORE if p in edges]
        if am:
            x = [edges[p] for p in common]; y = [am[p] for p in common]
            print("  %-10s against AlphaMissense Cliff delta   rho = %+.3f  over %d edges"
                  % (tag, spearman(x, y), len(common)))
        if hs:
            x = [edges[p] for p in common]; y = [-hs[p] for p in common]
            print("  %-10s against humsavar tolerance          rho = %+.3f  over %d edges"
                  % (tag, spearman(x, y), len(common)))

    # ------------------------------------------------------ the histidine routes
    print()
    print("=" * 100)
    print("the routes to histidine, which is what a pH-switch design installs")
    hdr = "%-9s %-12s %s" % ("route", "readout", "verdict")
    print(hdr)
    for tag, edges in results.items():
        for p in (("N", "H"), ("D", "H"), ("Y", "H")):
            if p in edges:
                print("  %s>%s     %-12s cliff %+.3f" % (p[0], p[1], tag, edges[p]))
    if am:
        for p in (("N", "H"), ("D", "H"), ("Y", "H")):
            if p in am:
                print("  %s>%s     %-12s cliff %+.3f" % (p[0], p[1], "alphamissense", am[p]))
    if hs:
        for p in (("N", "H"), ("D", "H")):
            if p in hs:
                print("  %s>%s     %-12s %.1f%% pathogenic" % (p[0], p[1], "humsavar", hs[p]))

    out = os.path.join(a.outdir, "ladder_ddg_core.csv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["edge", "direction"] + ["cliff_" + t for t, _ in sets]
                   + ["cliff_alphamissense_soluble", "humsavar_pathogenic_pct"])
        for p in CORE:
            w.writerow([p[0] + ">" + p[1],
                        "basic" if p in BASICDIR else "acidic" if (p[1], p[0]) in BASICDIR else "-"]
                       + ["%.3f" % results.get(t, {}).get(p, float("nan")) for t, _ in sets]
                       + ["%.3f" % am.get(p, float("nan")) if am else "",
                          "%.1f" % hs.get(p, float("nan")) if hs else ""])
    print("\nwrote", out)


if __name__ == "__main__":
    main()
