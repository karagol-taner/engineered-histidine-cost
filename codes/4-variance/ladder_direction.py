#!/usr/bin/env python3
import argparse, collections, csv, importlib.util, os, sys
import numpy as np


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


LM = _load(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        os.pardir, "ladder_megascale.py"), "ladder_megascale")


def _sets(targets, name):
    p = os.path.join(targets, name) if targets else ""
    if not p or not os.path.exists(p):
        p = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         os.pardir, "ladder_sets.py")
        name = "ladder_sets"
    return LM.load(p, name)


def bh(ps):
    m = len(ps)
    order = sorted(range(m), key=lambda i: ps[i])
    adj = [0.0] * m
    prev = 1.0
    for rk, i in enumerate(reversed(order), 1):
        v = min(prev, ps[i] * m / (m - rk + 1))
        adj[i] = prev = v
    return adj


def cliff_np(vals, base_sorted):
    nb = len(base_sorted)
    if nb == 0 or len(vals) == 0:
        return np.nan
    lt = np.searchsorted(base_sorted, vals, side="left").sum()
    gt = (nb * len(vals)) - np.searchsorted(base_sorted, vals, side="right").sum()
    return (gt - lt) / float(len(vals) * nb)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--megascale", default="../../data/megascale")
    ap.add_argument("--targets", default="")
    ap.add_argument("--source-null", default="../../data/ladder_source_null.csv")
    ap.add_argument("--out", default="../../data/ladder_direction.csv")
    ap.add_argument("--boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20261007)
    a = ap.parse_args()

    ls = _sets(a.targets, "ladder_edge_cost.py")
    ca = _sets(a.targets, "codon_audit.py")
    CORE = [tuple(p) for p in ls.CORE]
    BASIC = {tuple(p) for p in ls.BASICDIR}
    reach = ca.routes if hasattr(ca, "routes") else ls.routes
    if not callable(reach):
        sys.exit("FATAL: no reachability test in the substitution set module")

    raw = LM.read_megascale(a.megascale)
    if not raw:
        sys.exit("FATAL: no cDNA-proteolysis measurements read from %s" % a.megascale)
    flip = (sum(1 for r in raw if r[4] < 0) / len(raw)) > 0.5
    rows = [(d, a1, a2, (-v if flip else v)) for d, _, a1, a2, v, _ in raw]
    print("ladder_direction.py: each reciprocal pair, forward against reverse")
    print("cDNA-proteolysis single substitutions %d, costs mapped with higher = "
          "less stable" % len(rows))

    pairs = []
    for e in sorted(BASIC):
        rev = (e[1], e[0])
        if rev not in {tuple(p) for p in CORE}:
            sys.exit("FATAL: %s>%s has no reverse among the twelve" % e)
        pairs.append((e, rev))
    if len(pairs) != 6:
        sys.exit("FATAL: %d charge-raising directions, expected 6" % len(pairs))

    doms = sorted({d for d, _, _, _ in rows})
    by_dom = collections.defaultdict(lambda: collections.defaultdict(list))
    base_by_dom = collections.defaultdict(list)
    for d, a1, a2, v in rows:
        by_dom[d][(a1, a2)].append(v)
        if reach(a1, a2):
            base_by_dom[d].append(v)
    for d in doms:
        base_by_dom[d] = np.asarray(base_by_dom[d], dtype=float)
        for k in list(by_dom[d]):
            by_dom[d][k] = np.asarray(by_dom[d][k], dtype=float)
    print("%d domains, baseline rebuilt inside every resample" % len(doms))

    edges = sorted({e for p in pairs for e in p})

    def gaps(sel):
        base = np.sort(np.concatenate([base_by_dom[doms[i]] for i in sel]))
        cl = {}
        for e in edges:
            vs = [by_dom[doms[i]][e] for i in sel if e in by_dom[doms[i]]]
            cl[e] = cliff_np(np.concatenate(vs), base) if vs else np.nan
        return cl, [cl[f] - cl[r] for f, r in pairs]

    obs_cl, obs = gaps(range(len(doms)))
    rng = np.random.default_rng(a.seed)
    draws = rng.integers(0, len(doms), (a.boot, len(doms)))
    dist = np.empty((a.boot, len(pairs)))
    for b in range(a.boot):
        dist[b] = gaps(draws[b])[1]
        if (b + 1) % 1000 == 0:
            print("  %d of %d resamples" % (b + 1, a.boot), flush=True)

    ps, los, his = [], [], []
    for j in range(len(pairs)):
        d = dist[:, j]
        d = d[~np.isnan(d)]
        lo, hi = np.percentile(d, [2.5, 97.5])
        p = 2 * min((d <= 0).mean(), (d >= 0).mean())
        ps.append(max(p, 1.0 / a.boot))
        los.append(lo)
        his.append(hi)
    adj = bh(ps)

    src = {}
    if os.path.exists(a.source_null):
        for r in csv.DictReader(open(a.source_null, encoding="utf-8")):
            src[r["pair"]] = float(r["residual"])

    out = [("pair", "forward", "reverse", "gap", "lo", "hi", "p", "bh_p",
            "gap_after_source_null")]
    print("\n%-11s %8s %8s %8s %-20s %9s %9s %10s"
          % ("pair", "fwd", "rev", "gap", "95% interval", "p", "BH p", "gap, null"))
    for j, (f, r) in enumerate(pairs):
        fe, re_ = "%s>%s" % f, "%s>%s" % r
        g = obs[j]
        sg = (abs(src[fe] - src[re_]) if fe in src and re_ in src else float("nan"))
        print("%-11s %+8.3f %+8.3f %+8.3f  %+6.3f to %+-7.3f %9.4f %9.4f %10.3f"
              % (fe + " / " + re_, obs_cl[f], obs_cl[r], g, los[j], his[j],
                 ps[j], adj[j], sg))
        out.append((fe + " / " + re_, "%+.6f" % obs_cl[f], "%+.6f" % obs_cl[r],
                    "%+.6f" % g, "%+.6f" % los[j], "%+.6f" % his[j],
                    "%.6f" % ps[j], "%.6f" % adj[j],
                    "" if sg != sg else "%.6f" % sg))

    mean_all = sum(abs(x) for x in obs) / len(obs)
    hisx = [j for j, (f, _) in enumerate(pairs) if f[1] == "H"]
    rest = [j for j in range(len(pairs)) if j not in hisx]
    mh = sum(abs(obs[j]) for j in hisx) / len(hisx)
    mr = sum(abs(obs[j]) for j in rest) / len(rest)
    sn = [abs(src["%s>%s" % f] - src["%s>%s" % r]) for f, r in pairs
          if "%s>%s" % f in src and "%s>%s" % r in src]
    print("\nmean gap over the six pairs            %.3f" % mean_all)
    print("the %d pairs that install a histidine   %.3f" % (len(hisx), mh))
    print("the other %d pairs                      %.3f" % (len(rest), mr))
    if len(sn) == len(pairs):
        print("mean gap after the source-residue null %.3f, %.0f per cent of the measured"
              % (sum(sn) / len(sn), 100.0 * sum(sn) / sum(abs(x) for x in obs)))
    out.append(("mean, six pairs", "", "", "%+.6f" % mean_all, "", "", "", "", ""))
    out.append(("mean, histidine pairs", "", "", "%+.6f" % mh, "", "", "", "", ""))
    out.append(("mean, other pairs", "", "", "%+.6f" % mr, "", "", "", "", ""))
    if len(sn) == len(pairs):
        out.append(("mean after source null", "", "", "%+.6f" % (sum(sn) / len(sn)),
                    "", "", "", "", ""))

    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        csv.writer(fh).writerows(out)
    print("\nwrote %s" % a.out)


if __name__ == "__main__":
    main()
