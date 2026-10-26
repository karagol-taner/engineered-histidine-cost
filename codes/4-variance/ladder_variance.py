#!/usr/bin/env python3
"""Where the per-protein effect of a substitution lives.

Separates the variance of the per-protein effect into a between-substitution and
a protein-by-substitution component by restricted maximum likelihood, and
predicts held-out domains from the substitution, from the position, and from both.

    python3 ladder_variance.py --help
"""
import argparse, csv, glob, gzip, importlib.util, os, re, sys
import numpy as np

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


def rank(a):
    o = np.argsort(a, kind="mergesort")
    r = np.empty(len(a), dtype=np.float64)
    r[o] = np.arange(1, len(a) + 1)
    # average ties
    s = a[o]
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        if j > i:
            r[o[i:j + 1]] = (i + j) / 2.0 + 1
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = rank(np.asarray(x, float)), rank(np.asarray(y, float))
    return float(np.corrcoef(rx, ry)[0, 1])


def pearson(x, y):
    return float(np.corrcoef(np.asarray(x, float), np.asarray(y, float))[0, 1])


# ------------------------------------------------- variance decomposition
def _reml(d, v, gidx, ng, t_sub, t_int):
    """Restricted log-likelihood of the two-component model."""
    ld = 0.0
    a = b = 0.0            # for the GLS mean
    parts = []
    for g in range(ng):
        m = gidx == g
        if not m.any():
            continue
        D = t_int + v[m]
        w = 1.0 / D
        Sw = w.sum()
        f = 1.0 + t_sub * Sw
        ld += np.log(D).sum() + np.log(f)
        swd = (w * d[m]).sum()
        a += swd / f
        b += Sw / f
        parts.append((m, w, Sw, f))
    if b <= 0:
        return -np.inf
    mu = a / b
    quad = 0.0
    for m, w, Sw, f in parts:
        r = d[m] - mu
        swr = (w * r).sum()
        quad += (w * r * r).sum() - t_sub * swr * swr / f
    return -0.5 * (ld + quad) - 0.5 * np.log(b)


def decompose(cells, label, note=""):
    """cells: (group_key, substitution, effect, sampling_var, n)."""
    subs = sorted({c[1] for c in cells})
    if len(subs) < 3 or len(cells) < 20:
        print("\n--- %s ---\n    too few cells to decompose (%d cells, %d substitutions)"
              % (label, len(cells), len(subs)))
        return None
    d = np.array([c[2] for c in cells], float)
    v = np.array([c[3] for c in cells], float)
    g = np.array([subs.index(c[1]) for c in cells])
    ng = len(subs)

    spread = max(d.var(ddof=1), 1e-8)
    grid = np.concatenate([[0.0], np.exp(np.linspace(np.log(spread * 1e-4),
                                                     np.log(spread * 4), 48))])
    best = (-np.inf, 0.0, 0.0)
    for ts in grid:
        for ti in grid:
            ll = _reml(d, v, g, ng, ts, ti)
            if ll > best[0]:
                best = (ll, ts, ti)
    _, ts, ti = best
    for _ in range(4):                      # local refinement
        for scale in (0.5, 0.25):
            cand = [(max(0.0, ts * (1 + k * scale)), max(0.0, ti * (1 + j * scale)))
                    for k in (-1, 0, 1) for j in (-1, 0, 1)]
            for a2, b2 in cand:
                ll = _reml(d, v, g, ng, a2, b2)
                if ll > best[0]:
                    best = (ll, a2, b2)
            _, ts, ti = best
    tot = ts + ti
    share = ts / tot if tot > 0 else float("nan")
    print("\n--- %s ---" % label)
    if note:
        print("    %s" % note)
    print("    cells %d over %d substitutions and %d proteins, median cell n = %d"
          % (len(cells), ng, len({c[0] for c in cells}), int(np.median([c[4] for c in cells]))))
    print("    between-substitution variance    %.5f  (sd %.3f)" % (ts, ts ** 0.5))
    print("    protein-by-substitution variance %.5f  (sd %.3f)" % (ti, ti ** 0.5))
    print("    mean sampling variance           %.5f  (sd %.3f)" % (v.mean(), v.mean() ** 0.5))
    print("    SHARE ATTRIBUTABLE TO THE SUBSTITUTION: %.0f%%" % (100 * share))
    return dict(tau2_sub=ts, tau2_int=ti, samp=float(v.mean()), share=share,
                cells=len(cells), subs=ng, groups=len({c[0] for c in cells}))


def transfer(cells, label):
    """Leave-one-protein-out: predict each cell's effect from the SAME substitution
    measured in every OTHER protein."""
    if len(cells) < 50:
        print("\n--- %s ---\n    too few cells" % label); return None
    import collections
    bys = collections.defaultdict(list)
    for k, s, e, vv, n in cells:
        bys[s].append((k, e, vv, n))
    pred, obs = [], []
    for s, rows in bys.items():
        if len(rows) < 3:
            continue
        tot = sum(r[1] / max(r[2], 1e-9) for r in rows)
        wt = sum(1.0 / max(r[2], 1e-9) for r in rows)
        for k, e, vv, n in rows:
            w = 1.0 / max(vv, 1e-9)
            den = wt - w
            if den <= 0:
                continue
            pred.append((tot - e * w) / den)      # inverse-variance mean of the others
            obs.append(e)
    if len(pred) < 50:
        print("\n--- %s ---\n    too few usable cells" % label); return None
    r = pearson(pred, obs)
    print("\n--- %s ---" % label)
    print("    leave-one-protein-out, %d cells" % len(pred))
    print("    Pearson r = %+.3f, R squared = %.3f, Spearman rho = %+.3f"
          % (r, r * r, spearman(pred, obs)))
    return dict(r=r, r2=r * r, rho=spearman(pred, obs), n=len(pred))


# ---------------------------------------------------------------- readers
def read_megascale(d, core):
    rows = []
    for p in sorted(glob.glob(os.path.join(d, "*.csv"))):
        for r in csv.DictReader(open(p, newline="")):
            m = MUT.match((r.get("mut_type") or "").strip())
            if not m:
                continue
            a1, pos, a2 = m.groups()
            if a1 not in AA or a2 not in AA or a1 == a2:
                continue
            try:
                v = float(r.get("ddG_ML"))
            except (TypeError, ValueError):
                continue
            rows.append(((r.get("WT_name") or "").strip(), int(pos), a1, a2, v))
    if not rows:
        sys.exit("no megascale rows read from %s" % d)
    neg = sum(1 for r in rows if r[4] < 0) / len(rows)
    if neg > 0.5:
        rows = [(w, p, a, b, -v) for w, p, a, b, v in rows]
    pro = [v for _, _, _, b, v in rows if b == "P"]
    med = float(np.median([r[4] for r in rows]))
    if len(pro) < 30 or np.median(pro) <= med:
        sys.exit("FATAL: proline check failed on the megascale set")
    print("megascale: %d measurements, %.1f%% negative before mapping, proline median "
          "%+.3f against %+.3f overall" % (len(rows), 100 * neg, float(np.median(pro)), med))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="",
                    help="directory holding the method's own "
                         "modules; omit to use the bundled "
                         "codes/ladder_sets.py")
    ap.add_argument("--megascale", default="../../data/megascale")
    ap.add_argument("--cohort", nargs="*",
                    default=["../../data/alphamissense_transmembrane/*.csv"])
    ap.add_argument("--cohort-label", default="transmembrane")
    ap.add_argument("--min-cell", type=int, default=3)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=20261007)
    ap.add_argument("--out", default="../../data/ladder_variance.csv")
    a = ap.parse_args()

    lec = load(tpath(a.targets, "ladder_edge_cost.py"), "lec")
    ca = load(tpath(a.targets, "codon_audit.py"), "ca")
    CORE = lec.CORE
    CORE_S = {p[0] + ">" + p[1] for p in CORE}
    snv_ok = lambda x, y: bool(ca.routes(x, y))

    print("ladder_variance.py\n")
    print("=" * 100)
    print("1  VARIANCE DECOMPOSITION")

    rows = read_megascale(a.megascale, CORE_S)

    # per-domain baseline (single-nucleotide-reachable substitutions)
    base = {}
    for w, pos, a1, a2, v in rows:
        if snv_ok(a1, a2):
            base.setdefault(w, []).append(v)
    bmean = {k: float(np.mean(vs)) for k, vs in base.items() if len(vs) >= 50}

    def cells_from(rows, keep):
        acc = {}
        for w, pos, a1, a2, v in rows:
            s = a1 + ">" + a2
            if w not in bmean or (keep and s not in keep):
                continue
            acc.setdefault((w, s), []).append(v)
        out = []
        for (w, s), vs in acc.items():
            if len(vs) < a.min_cell:
                continue
            vs = np.array(vs)
            out.append((w, s, float(vs.mean()) - bmean[w],
                        float(vs.var(ddof=1)) / len(vs), len(vs)))
        return out

    res = {}
    mc = cells_from(rows, CORE_S)
    res["megascale_core"] = decompose(
        mc, "measured folding stability, the twelve core substitutions",
        "cell value is the mean cost of that substitution in that domain, minus the domain's own baseline")
    res["megascale_core_transfer"] = transfer(
        mc, "transferability, measured folding stability, the twelve core substitutions")
    ma = cells_from(rows, None)
    res["megascale_all"] = decompose(ma, "measured folding stability, all substitution classes")
    res["megascale_all_transfer"] = transfer(ma, "transferability, all substitution classes")

    # ---- the same decomposition on the variant cohort
    paths = []
    for p in a.cohort:
        paths.extend(sorted(glob.glob(os.path.expanduser(p))) or [os.path.expanduser(p)])
    acc, bacc = {}, {}
    for path in paths:
        op = gzip.open(path, "rt") if path.endswith(".gz") else open(path, newline="")
        with op as f:
            rd = csv.reader(f)
            header = [h.lstrip("# ").strip() for h in next(rd)]
            sc = next((c for c in lec.SCHEMAS if c["aa1"] in header and c["score"] in header), None)
            if sc is None:
                sys.exit("unrecognised header in %s" % path)
            iacc = next((header.index(k) for k in ("acc", "Uniprot ACC") if k in header), None)
            i1, i2 = header.index(sc["aa1"]), header.index(sc["aa2"])
            isc = header.index(sc["score"])
            isnv = header.index(sc["snv"]) if sc["snv"] in header else None
            for r in rd:
                if len(r) <= isc:
                    continue
                try:
                    v = float(r[isc])
                except ValueError:
                    continue
                w = r[iacc]
                if (isnv is None) or (r[isnv].strip().lower() == "y"):
                    bacc.setdefault(w, []).append(v)
                s = r[i1] + ">" + r[i2]
                if s in CORE_S:
                    acc.setdefault((w, s), []).append(v)
    bm = {k: float(np.mean(vs)) for k, vs in bacc.items() if len(vs) >= 200}
    cc = []
    for (w, s), vs in acc.items():
        if w not in bm or len(vs) < a.min_cell:
            continue
        vs = np.array(vs)
        cc.append((w, s, float(vs.mean()) - bm[w], float(vs.var(ddof=1)) / len(vs), len(vs)))
    res["cohort_core_transfer"] = transfer(
        cc, "transferability, %s cohort, the twelve core substitutions" % a.cohort_label)
    res["cohort_core"] = decompose(
        cc, "%s cohort, the twelve core substitutions" % a.cohort_label,
        "cell value is the mean pathogenicity score of that substitution in that protein, "
        "minus the protein's own baseline mean")

    # =================================================== 2 head to head
    print("\n" + "=" * 100)
    print("2  THE SUBSTITUTION AGAINST THE POSITION")
    print("Held-out domains, %d folds. The position signal is the mean measured cost of" % a.folds)
    print("every other substitution at the same position, which bounds what any")
    print("position-specific predictor could extract from that position.\n")

    dom = np.array([r[0] for r in rows])
    key = np.array([r[0] + "|" + str(r[1]) for r in rows])
    sub = np.array([r[2] + ">" + r[3] for r in rows])
    y = np.array([r[4] for r in rows], dtype=np.float64)

    uk, kinv = np.unique(key, return_inverse=True)
    psum = np.bincount(kinv, weights=y)
    pcnt = np.bincount(kinv)
    pos_loo = np.where(pcnt[kinv] > 1, (psum[kinv] - y) / np.maximum(pcnt[kinv] - 1, 1), np.nan)

    ud, dinv = np.unique(dom, return_inverse=True)
    rng = np.random.default_rng(a.seed)
    fold = rng.integers(0, a.folds, len(ud))[dinv]

    us, sinv = np.unique(sub, return_inverse=True)
    out_rows = []
    for name, mask in (("all substitution classes", np.ones(len(y), bool)),
                       ("the twelve core substitutions", np.isin(sub, list(CORE_S)))):
        ps, pp, pb, n_used = [], [], [], 0
        for f in range(a.folds):
            tr = (fold != f) & mask & np.isfinite(pos_loo)
            te = (fold == f) & mask & np.isfinite(pos_loo)
            if te.sum() < 50 or tr.sum() < 50:
                continue
            m = np.full(len(us), np.nan)
            cnt = np.bincount(sinv[tr], minlength=len(us))
            tot = np.bincount(sinv[tr], weights=y[tr], minlength=len(us))
            ok = cnt >= 5
            m[ok] = tot[ok] / cnt[ok]
            sp = m[sinv[te]]
            good = np.isfinite(sp)
            if good.sum() < 50:
                continue
            yt = y[te][good]; sp = sp[good]; pl = pos_loo[te][good]
            X = np.column_stack([np.ones(len(sp)), sp, pl])
            mtr = m[sinv[tr]]
            g2 = np.isfinite(mtr)
            Xtr = np.column_stack([np.ones(g2.sum()), mtr[g2], pos_loo[tr][g2]])
            beta, *_ = np.linalg.lstsq(Xtr, y[tr][g2], rcond=None)
            both = X @ beta
            ps.append((spearman(sp, yt), spearman(pl, yt), spearman(both, yt)))
            pp.append((pearson(sp, yt), pearson(pl, yt), pearson(both, yt)))
            n_used += len(yt)
        if not ps:
            continue
        ps = np.array(ps); pp = np.array(pp)
        print("--- %s, %d held-out measurements ---" % (name, n_used))
        print("    %-34s %-12s %-12s" % ("predictor", "Spearman", "R squared"))
        for j, lab in enumerate(("substitution identity alone",
                                 "position signal alone",
                                 "both together")):
            print("    %-34s %+.3f        %.3f" % (lab, ps[:, j].mean(), pp[:, j].mean() ** 2))
        gain_sub = pp[:, 2].mean() ** 2 - pp[:, 1].mean() ** 2
        gain_pos = pp[:, 2].mean() ** 2 - pp[:, 0].mean() ** 2
        print("    what the substitution adds on top of the position: %+.3f R squared" % gain_sub)
        print("    what the position adds on top of the substitution: %+.3f R squared" % gain_pos)
        out_rows.append((name, n_used, ps[:, 0].mean(), ps[:, 1].mean(), ps[:, 2].mean(),
                         pp[:, 0].mean() ** 2, pp[:, 1].mean() ** 2, pp[:, 2].mean() ** 2,
                         gain_sub, gain_pos))
        print()

    with open(os.path.expanduser(a.out), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["analysis", "key", "value"])
        for k, r in res.items():
            if not r:
                continue
            for f, val in r.items():
                w.writerow([k, f, "%.6f" % val if isinstance(val, float) else val])
        for r in out_rows:
            for f, v in zip(("n", "spearman_sub", "spearman_pos", "spearman_both",
                             "r2_sub", "r2_pos", "r2_both", "gain_sub", "gain_pos"), r[1:]):
                w.writerow(["headtohead:" + r[0], f, "%.4f" % v if isinstance(v, float) else v])
    print("wrote", os.path.expanduser(a.out))


if __name__ == "__main__":
    main()
