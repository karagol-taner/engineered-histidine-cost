#!/usr/bin/env python3
"""Designed against natural domains in the cDNA-proteolysis set.

Classifies each construct by its name, scores every substitution separately in
the two classes against its own baseline, and bootstraps the difference.

    python3 ladder_designed.py --help
"""
import argparse, bisect, csv, glob, importlib.util, os, re
import statistics as st, sys
import numpy as np

MUT = re.compile(r'^([A-Z])(-?\d+)([A-Z])$')
AA = set("ACDEFGHIKLMNPQRSTVWY")
DESIGN = re.compile(r'(^|_)(HHH|HHHH|EHEE|HEEH|EEHEE)_|TrROS|Hall|_rd\d|^rd\d|_de_|design', re.I)
NATURAL = re.compile(r'^(v2[A-Za-z0-9]*_)?[0-9][A-Za-z0-9]{3}(\.pdb)?$', re.I)


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


def classify(name):
    base = name.strip()
    if NATURAL.match(base.replace(".pdb", "")):
        return "natural"
    if DESIGN.search(base):
        return "designed"
    return None


def cliff_factory(base):
    base = sorted(base)
    nb = len(base)

    def cliff(xs):
        gt = lt = 0
        for x in xs:
            lt += bisect.bisect_left(base, x)
            gt += nb - bisect.bisect_right(base, x)
        return (gt - lt) / (len(xs) * nb)
    return cliff, nb


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="",
                    help="directory holding the method's own "
                         "modules; omit to use the bundled "
                         "codes/ladder_sets.py")
    ap.add_argument("--megascale", default="../../data/megascale")
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20261007)
    ap.add_argument("--out", default="../../data/ladder_designed.csv")
    a = ap.parse_args()

    lec = load(tpath(a.targets, "ladder_edge_cost.py"), "lec")
    ca = load(tpath(a.targets, "codon_audit.py"), "ca")
    CORE, EXT, BASICDIR = lec.CORE, lec.EXT, lec.BASICDIR
    snv_ok = lambda p: bool(ca.routes(p[0], p[1]))
    WANT = list(CORE) + [('Y', 'H'), ('H', 'Y'), ('K', 'R'), ('N', 'K')]

    print("ladder_designed.py\n")

    rows, unknown = [], {}
    for p in sorted(glob.glob(os.path.join(a.megascale, "*.csv"))):
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh):
                m = MUT.match((r.get("mut_type") or "").strip())
                if not m:
                    continue
                a1, _, a2 = m.groups()
                if a1 not in AA or a2 not in AA or a1 == a2:
                    continue
                try:
                    v = float(r.get("ddG_ML"))
                except (TypeError, ValueError):
                    continue
                w = (r.get("WT_name") or "").strip()
                c = classify(w)
                if c is None:
                    unknown[w] = unknown.get(w, 0) + 1
                    continue
                rows.append((a1, a2, v, c, w))
    if unknown:
        print("FATAL: %d constructs match neither the natural nor the designed "
              "pattern:" % len(unknown))
        for k, n in sorted(unknown.items())[:12]:
            print("   %-40s %d rows" % (k, n))
        sys.exit("refusing to classify by exclusion")
    if not rows:
        sys.exit("no rows read from %s" % a.megascale)

    neg = sum(1 for _, _, v, _, _ in rows if v < 0) / len(rows)
    flip = neg > 0.5
    rows = [(x, y, (-v if flip else v), c, w) for x, y, v, c, w in rows]
    print("sign: %.1f%% of values negative, so %s is destabilising; mapped to cost "
          "with higher = less stable" % (100 * neg, "negative" if flip else "positive"))
    med = st.median([v for _, _, v, _, _ in rows])
    pro = [v for _, y, v, _, _ in rows if y == 'P']
    if len(pro) < 30 or st.median(pro) <= med:
        sys.exit("FATAL: the proline check fails (n=%d), so the sign mapping is "
                 "not trustworthy." % len(pro))
    print("sign check, proline: median cost %+.3f overall against %+.3f for X>P "
          "(n=%d)\n" % (med, st.median(pro), len(pro)))

    doms = {}
    for _, _, _, c, w in rows:
        doms.setdefault(c, set()).add(w)
    for c in ("natural", "designed"):
        print("%-9s %6d measurements over %3d domains"
              % (c, sum(1 for r in rows if r[3] == c), len(doms.get(c, ()))))

    res, bootci = {}, {}
    rng = np.random.default_rng(a.seed)
    for c in ("natural", "designed"):
        sub = [(x, y, v, w) for x, y, v, cc, w in rows if cc == c]
        base = [v for x, y, v, _ in sub if snv_ok((x, y))]
        cliff, nb = cliff_factory(base)
        print("\n" + "=" * 96)
        print("%s domains: baseline is every single-nucleotide-reachable "
              "substitution, n=%d, mean %+.3f kcal/mol, median %+.3f"
              % (c.upper(), nb, st.mean(base), st.median(base)))
        print("%-7s %-7s %-8s %-10s %-10s %-9s %s"
              % ("edge", "dir", "n", "mean", "median", "cliff", "95% interval"))
        per = {}
        for x, y, v, w in sub:
            per.setdefault((x, y), []).append((v, w))
        dom_list = sorted(doms[c])
        dom_ix = {d: i for i, d in enumerate(dom_list)}
        nd = len(dom_list)
        idx = rng.integers(0, nd, (a.boot, nd))
        out = []
        for p in WANT:
            vv = per.get(p)
            if not vv or len(vv) < 10:
                continue
            vals = [v for v, _ in vv]
            d = cliff(vals)
            by = {}
            for v, w in vv:
                by.setdefault(dom_ix[w], []).append(v)
            bb = []
            for b in range(a.boot):
                xs = [v for j in idx[b] if j in by for v in by[j]]
                if xs:
                    bb.append(cliff(xs))
            lo, hi = (np.percentile(bb, 2.5), np.percentile(bb, 97.5)) if bb else (float('nan'),) * 2
            res.setdefault(c, {})[p] = (len(vals), st.mean(vals), st.median(vals), d)
            bootci.setdefault(c, {})[p] = (lo, hi)
            out.append((-d, "%-7s %-7s %-8d %+-10.3f %+-10.3f %+-9.3f %+.3f to %+.3f"
                        % (p[0] + ">" + p[1],
                           "basic" if p in BASICDIR else
                           "acidic" if (p[1], p[0]) in BASICDIR else "-",
                           len(vals), st.mean(vals), st.median(vals), d, lo, hi)))
        for _, l in sorted(out):
            print(l)

    print("\n" + "=" * 96)
    print("DESIGNED AGAINST NATURAL, the comparison this was built for")
    common = [p for p in CORE if p in res.get("natural", {}) and p in res.get("designed", {})]
    rho = spearman([res["natural"][p][3] for p in common],
                   [res["designed"][p][3] for p in common])
    print("Spearman rho between the two orderings: %+.3f over %d core substitutions"
          % (rho, len(common)))
    print("\n%-7s %-24s %-24s %s" % ("edge", "natural cliff (n)", "designed cliff (n)", "difference"))
    for p in sorted(common, key=lambda q: -res["natural"][q][3]):
        n_, d_ = res["natural"][p], res["designed"][p]
        print("%-7s %+-8.3f (n=%-6d)   %+-8.3f (n=%-6d)   %+.3f"
              % (p[0] + ">" + p[1], n_[3], n_[0], d_[3], d_[0], d_[3] - n_[3]))
    print("\nthe routes to histidine, in kcal/mol and on the Cliff scale")
    for p in (("N", "H"), ("D", "H"), ("Y", "H")):
        if p in res.get("natural", {}) and p in res.get("designed", {}):
            n_, d_ = res["natural"][p], res["designed"][p]
            ln, hn = bootci["natural"][p]
            ld, hd = bootci["designed"][p]
            print("  %s>%s  natural  %+.3f kcal/mol, cliff %+.3f (%+.3f to %+.3f), n=%d"
                  % (p[0], p[1], n_[1], n_[3], ln, hn, n_[0]))
            print("       designed %+.3f kcal/mol, cliff %+.3f (%+.3f to %+.3f), n=%d"
                  % (d_[1], d_[3], ld, hd, d_[0]))

    with open(os.path.expanduser(a.out), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["edge", "direction", "class", "n", "mean_kcal", "median_kcal",
                    "cliff", "boot_lo", "boot_hi"])
        for c in ("natural", "designed"):
            for p in WANT:
                if p not in res.get(c, {}):
                    continue
                n_, m_, md_, d_ = res[c][p]
                lo, hi = bootci[c][p]
                w.writerow([p[0] + ">" + p[1],
                            "basic" if p in BASICDIR else
                            "acidic" if (p[1], p[0]) in BASICDIR else "-",
                            c, n_, "%.3f" % m_, "%.3f" % md_, "%.3f" % d_,
                            "%.3f" % lo, "%.3f" % hi])
    print("\nwrote", os.path.expanduser(a.out))


if __name__ == "__main__":
    main()
