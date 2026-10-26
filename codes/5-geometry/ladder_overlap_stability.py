#!/usr/bin/env python3
"""Shape overlap against measured folding cost, over every substitution class.

Scores every substitution class in the cDNA-proteolysis set against the
whole-set baseline and correlates it with the overlap of its residue pair, with
partial correlations holding size constant.

    python3 ladder_overlap_stability.py --help
"""
import argparse, bisect, csv, glob, math, os, random, re, statistics as st, sys

MUT = re.compile(r"^([A-Z])(\d+)([A-Z])$")
AA = set("ACDEFGHIKLMNPQRSTVWY")
TWELVE = ["Q>E", "K>Q", "E>Q", "Q>K", "N>D", "N>H",
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


def perm_p(x, y, reps, seed):
    r0 = abs(spear(x, y))
    rnd = random.Random(seed)
    y2 = list(y)
    c = 0
    for _ in range(reps):
        rnd.shuffle(y2)
        if abs(spear(x, y2)) >= r0 - 1e-12:
            c += 1
    return (c + 1) / (reps + 1)


def partial(x, y, z):
    rx, ry, rz = rank(x), rank(y), rank(z)
    a, b, c = pear(rx, ry), pear(rx, rz), pear(ry, rz)
    return (a - b * c) / math.sqrt((1 - b ** 2) * (1 - c ** 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--overlaps", default="../../data/ladder_overlap_all_pairs.csv")
    ap.add_argument("--megascale", default="../../data/megascale")
    ap.add_argument("--min-n", type=int, default=3)
    ap.add_argument("--drop", default="C",
                    help="residue types with too few observed side chains to bin")
    ap.add_argument("--perm", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20261007)
    ap.add_argument("--classes-out", default="../../data/ladder_class_costs.csv")
    ap.add_argument("--out", default="../../data/ladder_overlap_stability.csv")
    a = ap.parse_args()
    print("ladder_overlap_stability.py\n")

    ov, dre, dat = {}, {}, {}
    nres = {}
    for r in csv.DictReader(open(a.overlaps, newline="")):
        ov[r["pair"]] = float(r["overlap"])
        dre[r["pair"]] = float(r["d_reach"])
        dat[r["pair"]] = float(r["d_atoms"])
        nres[r["pair"][0]] = int(r["n_src"])
    if not ov:
        die("no overlaps in %s" % a.overlaps)

    vals, allv = {}, []
    for f in sorted(glob.glob(os.path.join(a.megascale, "*.csv"))):
        for r in csv.DictReader(open(f, newline="")):
            m = MUT.match((r.get("mut_type") or "").strip())
            if not m:
                continue
            w, t = m.group(1), m.group(3)
            if w not in AA or t not in AA or w == t:
                continue
            try:
                v = -float(r["ddG_ML"])
            except (TypeError, ValueError):
                continue
            vals.setdefault(w + ">" + t, []).append(v)
            allv.append(v)
    if not allv:
        die("no measurements read from %s" % a.megascale)
    allv.sort()
    nb = len(allv)

    def cliff(xs):
        gt = lt = 0
        for x in xs:
            lt += bisect.bisect_left(allv, x)
            gt += nb - bisect.bisect_right(allv, x)
        return (gt - lt) / (len(xs) * nb)

    rows = [("pair", "n", "mean_cost", "median_cost", "cliff")]
    cl = {}
    for k in sorted(vals):
        if len(vals[k]) < a.min_n:
            continue
        cl[k] = cliff(vals[k])
        rows.append((k, len(vals[k]), "%.4f" % st.mean(vals[k]),
                     "%.4f" % st.median(vals[k]), "%+.4f" % cl[k]))
    with open(a.classes_out, "w", newline="") as fh:
        csv.writer(fh).writerows(rows)
    print("%d measurements, %d classes at %d or more" % (nb, len(cl), a.min_n))

    drop = set(a.drop)
    pairs = [p for p in ov if p in cl and p[0] not in drop and p[2] not in drop]
    und = {}
    for p in pairs:
        und.setdefault("".join(sorted([p[0], p[2]])), []).append(cl[p])
    uk = [k for k in und if len(und[k]) == 2]
    uo = [ov[k[0] + ">" + k[1]] for k in uk]
    uc = [sum(und[k]) / 2 for k in uk]
    ua = [abs(dat[k[0] + ">" + k[1]]) for k in uk]
    ur = [abs(dre[k[0] + ">" + k[1]]) for k in uk]
    o = [ov[p] for p in pairs]
    c = [cl[p] for p in pairs]
    tw = [p for p in TWELVE if p in ov]

    res = [("quantity", "value", "n", "p")]
    res.append(("overlap vs Cliff delta, ordered classes", "%+.3f" % spear(o, c),
                len(pairs), "%.5f" % perm_p(o, c, a.perm, a.seed)))
    res.append(("overlap vs mean Cliff delta, residue pairs", "%+.3f" % spear(uo, uc),
                len(uk), "%.5f" % perm_p(uo, uc, a.perm, a.seed)))
    res.append(("the same, change in heavy atoms held constant",
                "%+.3f" % partial(uo, uc, ua), len(uk), ""))
    res.append(("the same, change in reach held constant",
                "%+.3f" % partial(uo, uc, ur), len(uk), ""))
    res.append(("|change in heavy atoms| vs mean Cliff delta", "%+.3f" % spear(ua, uc),
                len(uk), ""))
    res.append(("|change in reach| vs mean Cliff delta", "%+.3f" % spear(ur, uc),
                len(uk), ""))
    res.append(("overlap range, all classes", "%.3f to %.3f" % (min(o), max(o)),
                len(pairs), ""))
    res.append(("overlap range, the twelve",
                "%.3f to %.3f" % (min(ov[p] for p in tw), max(ov[p] for p in tw)), 12, ""))
    res.append(("residue types measured", str(len(nres)), str(sum(nres.values())), ""))
    res.append(("residue types used in the test",
                str(len(set(p[0] for p in pairs))),
                str(sum(v for k, v in nres.items() if k not in drop)), ""))
    with open(a.out, "w", newline="") as fh:
        csv.writer(fh).writerows(res)
    w = max(len(r[0]) for r in res)
    for r in res[1:]:
        print("  %-*s %-16s %-6s %s" % (w, r[0], r[1], r[2], r[3]))
    print("\nwrote %s and %s" % (a.classes_out, a.out))


if __name__ == "__main__":
    main()
