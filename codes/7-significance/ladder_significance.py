#!/usr/bin/env python3
"""Significance of each substitution and each rank agreement, with correction.

Tests every substitution against its baseline by cluster bootstrap over
proteins and every rank agreement by permutation, with Benjamini-Hochberg
applied across the twelve and across the eight comparisons.

    python3 ladder_significance.py
"""
import csv, os, random, sys
import numpy as np

if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
    print(__doc__.strip())
    sys.exit(0)

HERE = os.path.dirname(os.path.abspath(__file__))


def _root():
    d = HERE
    for _ in range(4):
        if os.path.isdir(os.path.join(d, "data")):
            return d
        d = os.path.dirname(d)
    sys.exit("FATAL: cannot find the data directory from %s" % HERE)


D = _root()
BOOT, PERM, SEED = 10000, 200000, 20261007
CORE = ["Q>E", "K>Q", "E>Q", "Q>K", "N>D", "N>H",
        "H>N", "K>E", "D>N", "E>K", "H>D", "D>H"]


def die(m):
    sys.exit("FATAL: " + m)


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


def pear(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    den = (sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y)) ** 0.5
    return num / den if den else float("nan")


def spear(x, y):
    return pear(rank(x), rank(y))


def perm_p(x, y, seed=SEED):
    """Two-sided permutation p for a Spearman rho, exact when n! is small."""
    import itertools
    n, r0 = len(x), abs(spear(x, y))
    fact = 1
    for i in range(2, n + 1):
        fact *= i
    if fact <= 500000:
        c = sum(1 for p in itertools.permutations(y)
                if abs(spear(x, list(p))) >= r0 - 1e-12)
        return c / fact, "exact over %d permutations" % fact
    rnd = random.Random(seed)
    y2 = list(y)
    c = 0
    for _ in range(PERM):
        rnd.shuffle(y2)
        if abs(spear(x, y2)) >= r0 - 1e-12:
            c += 1
    return (c + 1) / (PERM + 1), "%d random permutations" % PERM


def bh(ps):
    """Benjamini-Hochberg adjusted p-values, in the order given."""
    m = len(ps)
    order = sorted(range(m), key=lambda i: ps[i])
    adj = [0.0] * m
    prev = 1.0
    for rk, i in enumerate(reversed(order), 1):
        v = min(prev, ps[i] * m / (m - rk + 1))
        adj[i] = prev = v
    return adj


# ------------------------------------------- per-cell cluster bootstrap
def cell_bootstrap(path, label, out_rows):
    cells = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            cells.setdefault(r["protein_index"], {})[r["edge"]] = (
                int(r["n"]), float(r["cliff_within_protein"]))
    prots = sorted(cells)
    rng = np.random.default_rng(SEED)
    draws = rng.integers(0, len(prots), (BOOT, len(prots)))

    def pooled(sel, e):
        num = den = 0.0
        for i in sel:
            d = cells[prots[i]].get(e)
            if d:
                num += d[0] * d[1]
                den += d[0]
        return num / den if den else np.nan

    base = list(range(len(prots)))
    print("\n%s: %d proteins, cluster bootstrap over proteins, %d resamples"
          % (label, len(prots), BOOT))
    print("%-6s %9s %-20s %10s %10s" % ("edge", "cliff", "95% interval", "p", "BH p"))
    obs, dist, ps = {}, {}, []
    for e in CORE:
        obs[e] = pooled(base, e)
        d = np.array([pooled(draws[b], e) for b in range(BOOT)])
        d = d[~np.isnan(d)]
        dist[e] = d
        p = 2 * min((d <= 0).mean(), (d >= 0).mean())
        ps.append(max(p, 1.0 / BOOT))
    adj = bh(ps)
    for e, p, a in zip(CORE, ps, adj):
        lo, hi = np.percentile(dist[e], [2.5, 97.5])
        print("%-6s %+9.3f  %+6.3f to %+-7.3f %10.4f %10.4f" % (e, obs[e], lo, hi, p, a))
        out_rows.append((label, e, "%+.3f" % obs[e], "%+.3f" % lo, "%+.3f" % hi,
                         "%.4f" % p, "%.4f" % a))
    d = dist["N>H"] - dist["D>H"]
    lo, hi = np.percentile(d, [2.5, 97.5])
    p = 2 * min((d <= 0).mean(), (d >= 0).mean())
    print("  N>H minus D>H: %+.3f, 95%% %+.3f to %+.3f, p = %.4f"
          % (obs["N>H"] - obs["D>H"], lo, hi, max(p, 1.0 / BOOT)))
    out_rows.append((label, "N>H-D>H", "%+.3f" % (obs["N>H"] - obs["D>H"]),
                     "%+.3f" % lo, "%+.3f" % hi, "%.4f" % max(p, 1.0 / BOOT), ""))
    return obs


rows = [("readout", "edge", "value", "lo", "hi", "p", "bh_p")]
print("=" * 78)
print("1  EACH SUBSTITUTION AGAINST ITS BASELINE")
for f, lab in (("ladder_percell_soluble.csv", "AlphaMissense soluble"),
               ("ladder_percell_transmembrane.csv", "AlphaMissense transmembrane")):
    p = os.path.join(D, "data", f)
    if not os.path.exists(p):
        die("missing %s" % f)
    cell_bootstrap(p, lab, rows)

print("\n" + "=" * 78)
print("2  RANK AGREEMENT BETWEEN READOUTS")
core = {}
with open(os.path.join(D, "data", "ladder_core_fullprec.csv"), newline="") as fh:
    for r in csv.DictReader(fh):
        core[r["edge"]] = r
TMcliff = {}
with open(os.path.join(D, "data", "ladder_robustness_transmembrane.csv"), newline="") as fh:
    for r in csv.DictReader(fh):
        TMcliff[r["edge"]] = float(r["pooled_cliff"])
vals = [float(core[e]["cliff_megascale"]) for e in CORE]
if len(set(vals)) != len(vals):
    die("the stability column holds a tie; the full-precision file is required here")


def col(name, neg=False):
    out = {}
    for e in CORE:
        v = core[e][name]
        if v in ("", "nan"):
            continue
        out[e] = -float(v) if neg else float(v)
    return out


series = {
    "megascale": col("cliff_megascale"),
    "fireprot": col("cliff_fireprot"),
    "s669": col("cliff_s669"),
    "AM soluble": col("cliff_alphamissense_soluble"),
    "AM transmembrane": {e: TMcliff[e] for e in CORE if e in TMcliff},
    "humsavar tolerance": col("humsavar_pathogenic_pct", neg=True),
}
pairs = [("megascale", "AM soluble"), ("megascale", "humsavar tolerance"),
         ("fireprot", "AM soluble"), ("fireprot", "humsavar tolerance"),
         ("s669", "AM soluble"), ("s669", "humsavar tolerance"),
         ("humsavar tolerance", "AM soluble"), ("AM soluble", "AM transmembrane"),
         ("megascale", "AM transmembrane"), ("humsavar tolerance", "AM transmembrane"),
         ("megascale", "fireprot")]
print("%-38s %7s %5s %10s  %s" % ("comparison", "rho", "n", "p", "basis"))
pr = []
for a, b in pairs:
    sh = [e for e in CORE if e in series[a] and e in series[b]]
    if len(sh) < 4:
        continue
    x = [series[a][e] for e in sh]
    y = [series[b][e] for e in sh]
    p, how = perm_p(x, y)
    pr.append((a, b, spear(x, y), len(sh), p, how))
adj = bh([r[4] for r in pr])
for (a, b, rho, n, p, how), q in zip(pr, adj):
    print("%-38s %+7.3f %5d %10.5f  %s" % (a + " vs " + b, rho, n, p, how))
    rows.append(("rank agreement", a + " vs " + b, "%+.3f" % rho, "", "",
                 "%.5f" % p, "%.4f" % q))

print("\n" + "=" * 78)
print("3  GEOMETRY AGAINST THE ORDERING")
G = {}
with open(os.path.join(D, "data", "ladder_geometry.csv"), newline="") as fh:
    for r in csv.DictReader(fh):
        G.setdefault(r["analysis"], {})[r["key"]] = r["value"]
ms = {e: float(core[e]["cliff_megascale"]) for e in CORE}
x12 = [abs(float(G["d_reach"][e])) for e in CORE]
y12 = [ms[e] for e in CORE]
p12, how12 = perm_p(x12, y12)
print("  |change in reach| vs megascale Cliff, over the twelve")
print("    rho %+.3f, p = %.5f, %s" % (spear(x12, y12), p12, how12))
seen, xs, ys = set(), [], []
for e in CORE:
    k = "".join(sorted([e[0], e[2]]))
    if k in seen:
        continue
    seen.add(k)
    xs.append(abs(float(G["d_reach"][e])))
    ys.append((ms[e] + ms[e[2] + ">" + e[0]]) / 2.0)
pp, howp = perm_p(xs, ys)
print("  the same at pair level, where the descriptor actually varies")
print("    rho %+.3f over %d pairs, p = %.5f, %s" % (spear(xs, ys), len(xs), pp, howp))
rows.append(("geometry", "|d reach| vs cliff, twelve", "%+.3f" % spear(x12, y12),
             "", "", "%.5f" % p12, ""))
rows.append(("geometry", "|d reach| vs cliff, six pairs", "%+.3f" % spear(xs, ys),
             "", "", "%.5f" % pp, ""))

out = os.path.join(D, "data", "ladder_significance.csv")
with open(out, "w", newline="") as fh:
    csv.writer(fh).writerows(rows)
print("\nwrote %s" % out)
