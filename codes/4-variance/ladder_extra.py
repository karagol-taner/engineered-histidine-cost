#!/usr/bin/env python3
"""Four checks on the ordering.

Regresses the Cliff delta on charge, volume and hydropathy; recomputes it against
alternative baselines; compares topology classes; and counts transition against
transversion routes in observed human variation.

    python3 ladder_extra.py --help
"""
import argparse, bisect, csv, glob, importlib.util, os, re, sys
import statistics as st
import numpy as np

MUT = re.compile(r'^([A-Z])(-?\d+)([A-Z])$')
AA = set("ACDEFGHIKLMNPQRSTVWY")

# Zamyatnin residue volumes (cubic angstroms) and Kyte-Doolittle hydropathy.
VOL = dict(A=88.6, R=173.4, N=114.1, D=111.1, C=108.5, Q=143.8, E=138.4, G=60.1,
           H=153.2, I=166.7, L=166.7, K=168.6, M=162.9, F=189.9, P=112.7, S=89.0,
           T=116.1, W=227.8, Y=193.6, V=140.0)
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2,
          I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9,
          Y=-1.3, V=4.2)
DESIGN_ALPHA = re.compile(r'(^|_)(HHH|HHHH)_')
DESIGN_MIXED = re.compile(r'(^|_)(EHEE|HEEH|EEHEE)_')


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


def cliff_factory(base):
    base = sorted(base); nb = len(base)
    def cliff(xs):
        gt = lt = 0
        for x in xs:
            lt += bisect.bisect_left(base, x)
            gt += nb - bisect.bisect_right(base, x)
        return (gt - lt) / (len(xs) * nb)
    return cliff, nb


def spearman(x, y):
    def rk(v):
        s = sorted(range(len(v)), key=lambda i: v[i]); r=[0.0]*len(v); i=0
        while i < len(s):
            j=i
            while j+1 < len(s) and v[s[j+1]] == v[s[i]]: j+=1
            for k in range(i,j+1): r[s[k]] = (i+j)/2.0+1
            i=j+1
        return r
    rx, ry = rk(x), rk(y); n=len(x)
    mx, my = sum(rx)/n, sum(ry)/n
    num = sum((a-mx)*(b-my) for a,b in zip(rx,ry))
    den = (sum((a-mx)**2 for a in rx)*sum((b-my)**2 for b in ry))**0.5
    return num/den if den else float("nan")


def r2(y, X):
    X = np.column_stack([np.ones(len(y))] + list(X))
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ beta
    ss = ((y - y.mean()) ** 2).sum()
    return 1 - ((y - pred) ** 2).sum() / ss if ss > 0 else float("nan")


def read_megascale(d):
    rows = []
    for p in sorted(glob.glob(os.path.join(d, "*.csv"))):
        for r in csv.DictReader(open(p, newline="")):
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
            rows.append(((r.get("WT_name") or "").strip(), a1, a2, v))
    neg = sum(1 for r in rows if r[3] < 0) / len(rows)
    if neg > 0.5:
        rows = [(w, a, b, -v) for w, a, b, v in rows]
    pro = [v for _, _, b, v in rows if b == "P"]
    med = st.median([r[3] for r in rows])
    if len(pro) < 30 or st.median(pro) <= med:
        sys.exit("FATAL: proline check failed")
    print("megascale: %d measurements, proline median %+.3f against %+.3f overall"
          % (len(rows), st.median(pro), med))
    return rows


def parse_humsavar(path):
    """edge -> (LP/P, LB/B, US) from a ladder_humsavar run."""
    out, inblock = {}, False
    for line in open(path):
        if line.startswith("--- 12 core edges"):
            inblock = True; continue
        if inblock and line.startswith("---"):
            break
        f = line.split()
        if inblock and len(f) > 4 and ">" in f[0]:
            try:
                out[f[0]] = (int(f[2]), int(f[3]), int(f[4]))
            except ValueError:
                pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="",
                    help="directory holding the method's own "
                         "modules; omit to use the bundled "
                         "codes/ladder_sets.py")
    ap.add_argument("--megascale", default="../../data/megascale")
    ap.add_argument("--fireprot", default="../../data/fireprot")
    ap.add_argument("--humsavar-out", default="../../data/source_runs/ladder_humsavar.out")
    ap.add_argument("--min-n", type=int, default=3)
    ap.add_argument("--out", default="../../data/ladder_extra.csv")
    a = ap.parse_args()

    lec = load(tpath(a.targets, "ladder_edge_cost.py"), "lec")
    ca = load(tpath(a.targets, "codon_audit.py"), "ca")
    CORE, CHARGE = lec.CORE, lec.CHARGE
    CORE_S = [p[0] + ">" + p[1] for p in CORE]
    snv_ok = lambda x, y: bool(ca.routes(x, y))
    rows = read_megascale(a.megascale)
    out_rows = []

    # ---------------------------------------------------------------- 1 & 2
    per = {}
    for w, a1, a2, v in rows:
        per.setdefault(a1 + ">" + a2, []).append(v)
    base_snv = [v for w, a1, a2, v in rows if snv_ok(a1, a2)]
    base_all = [v for w, a1, a2, v in rows]
    c_snv, n_snv = cliff_factory(base_snv)
    c_all, n_all = cliff_factory(base_all)

    print("\n" + "=" * 96)
    print("1  IS THE ORDERING SIDE-CHAIN CHEMISTRY?")
    print("charge magnitude, side-chain volume and hydropathy, as magnitudes\n")
    for name, keys in (("the twelve core substitutions", CORE_S),
                       ("all substitution classes with at least %d measurements" % a.min_n,
                        [k for k, v in per.items() if len(v) >= a.min_n])):
        ks = [k for k in keys if k in per and len(per[k]) >= a.min_n]
        y = np.array([c_snv(per[k]) for k in ks])
        x1, x2 = [k.split(">")[0] for k in ks], [k.split(">")[1] for k in ks]
        dq = np.array([abs(CHARGE.get(b, 0) - CHARGE.get(x, 0)) if b in CHARGE and x in CHARGE
                       else abs((1 if b in "KR" else -1 if b in "DE" else 0.1 if b == "H" else 0) -
                                (1 if x in "KR" else -1 if x in "DE" else 0.1 if x == "H" else 0))
                       for x, b in zip(x1, x2)])
        dv = np.array([abs(VOL[b] - VOL[x]) for x, b in zip(x1, x2)])
        dh = np.array([abs(KD[b] - KD[x]) for x, b in zip(x1, x2)])
        print("--- %s, n = %d ---" % (name, len(ks)))
        print("    charge magnitude alone                  R squared %.3f" % r2(y, [dq]))
        print("    charge and volume                       R squared %.3f" % r2(y, [dq, dv]))
        print("    charge, volume and hydropathy           R squared %.3f" % r2(y, [dq, dv, dh]))
        print("    volume and hydropathy, no charge        R squared %.3f\n" % r2(y, [dv, dh]))
        out_rows.append(("descriptors:" + name, "n", len(ks)))
        out_rows.append(("descriptors:" + name, "r2_charge", "%.3f" % r2(y, [dq])))
        out_rows.append(("descriptors:" + name, "r2_charge_vol_hyd", "%.3f" % r2(y, [dq, dv, dh])))

    print("=" * 96)
    print("2  DOES THE ORDERING DEPEND ON THE BASELINE?\n")
    base_q = {}
    for w, a1, a2, v in rows:
        q = abs((1 if a2 in "KR" else -1 if a2 in "DE" else 0.1 if a2 == "H" else 0) -
                (1 if a1 in "KR" else -1 if a1 in "DE" else 0.1 if a1 == "H" else 0))
        base_q.setdefault(round(q, 1), []).append(v)
    def cliff_matched(k):
        x, b = k.split(">")
        q = abs((1 if b in "KR" else -1 if b in "DE" else 0.1 if b == "H" else 0) -
                (1 if x in "KR" else -1 if x in "DE" else 0.1 if x == "H" else 0))
        pool = base_q.get(round(q, 1), [])
        if len(pool) < 200:
            return float("nan")
        c, _ = cliff_factory(pool)
        return c(per[k])
    cols = {}
    for k in CORE_S:
        cols[k] = (c_snv(per[k]), c_all(per[k]), cliff_matched(k))
    print("%-7s %-12s %-12s %-12s" % ("edge", "SNV baseline", "all subs", "charge-matched"))
    for k in sorted(CORE_S, key=lambda e: -cols[e][0]):
        print("%-7s %+-12.3f %+-12.3f %s" % (k, cols[k][0], cols[k][1],
              ("%+.3f" % cols[k][2]) if cols[k][2] == cols[k][2] else "   .  "))
    ks = [k for k in CORE_S if cols[k][2] == cols[k][2]]
    r_all = spearman([cols[k][0] for k in CORE_S], [cols[k][1] for k in CORE_S])
    r_q = spearman([cols[k][0] for k in ks], [cols[k][2] for k in ks])
    print("\nrank agreement with the published baseline: all substitutions rho %+.3f, "
          "charge-matched rho %+.3f over %d edges" % (r_all, r_q, len(ks)))
    print("baselines: single-nucleotide n=%d, all n=%d" % (n_snv, n_all))
    out_rows += [("baseline", "rho_all_subs", "%.3f" % r_all),
                 ("baseline", "rho_charge_matched", "%.3f" % r_q)]

    # ------------------------------------------------------------------ 3
    print("\n" + "=" * 96)
    print("3  DOES THE ORDERING DEPEND ON THE FOLD?\n")
    cls = {}
    for w, a1, a2, v in rows:
        c = "all-alpha" if DESIGN_ALPHA.search(w) else "mixed alpha/beta" if DESIGN_MIXED.search(w) else None
        if c:
            cls.setdefault(c, []).append((a1, a2, v))
    for c in ("all-alpha", "mixed alpha/beta"):
        sub = cls.get(c, [])
        print("%-18s %6d measurements over %3d constructs"
              % (c, len(sub), len({w for w, a1, a2, v in rows
                                   if (DESIGN_ALPHA.search(w) if c == "all-alpha" else DESIGN_MIXED.search(w))})))
    got = {}
    for c in ("all-alpha", "mixed alpha/beta"):
        sub = cls.get(c, [])
        if not sub:
            continue
        b = [v for x, y, v in sub if snv_ok(x, y)]
        cf, nb = cliff_factory(b)
        pc = {}
        for x, y, v in sub:
            pc.setdefault(x + ">" + y, []).append(v)
        got[c] = {k: (cf(pc[k]), len(pc[k])) for k in CORE_S if k in pc and len(pc[k]) >= 10}
        print("   %s baseline n=%d, mean %+.3f" % (c, nb, st.mean(b)))
    common = [k for k in CORE_S if k in got.get("all-alpha", {}) and k in got.get("mixed alpha/beta", {})]
    if len(common) >= 4:
        print("\n%-7s %-22s %-22s" % ("edge", "all-alpha", "mixed alpha/beta"))
        for k in sorted(common, key=lambda e: -got["all-alpha"][e][0]):
            print("%-7s %+-8.3f (n=%-5d)      %+-8.3f (n=%-5d)"
                  % (k, got["all-alpha"][k][0], got["all-alpha"][k][1],
                     got["mixed alpha/beta"][k][0], got["mixed alpha/beta"][k][1]))
        rho = spearman([got["all-alpha"][k][0] for k in common],
                       [got["mixed alpha/beta"][k][0] for k in common])
        print("\nrank agreement between the two topology classes: rho %+.3f over %d edges"
              % (rho, len(common)))
        out_rows.append(("topology", "rho_alpha_vs_mixed", "%.3f" % rho))
    else:
        print("\nonly %d substitutions clear the threshold in both classes, too few to rank"
              % len(common))

    # FireProt secondary structure
    seen, fp = set(), []
    for p in sorted(glob.glob(os.path.join(a.fireprot, "*.csv"))):
        for r in csv.DictReader(open(p, newline="")):
            a1, a2 = r.get("wild_type"), r.get("mutation")
            try:
                v = float(r.get("ddG"))
            except (TypeError, ValueError):
                continue
            if not a1 or not a2 or a1 == a2:
                continue
            k = (r.get("uniprot_id"), r.get("position"), a1, a2, round(v, 3))
            if k in seen:
                continue
            seen.add(k)
            fp.append((a1, a2, v, (r.get("secondary_structure") or "").strip()))
    if fp:
        neg = sum(1 for _, _, v, _ in fp if v < 0) / len(fp)
        if neg > 0.5:
            fp = [(x, y, -v, s) for x, y, v, s in fp]
        from collections import Counter
        print("\nFireProtDB secondary structure of the substituted position: %s"
              % dict(Counter(s for _, _, _, s in fp).most_common(6)))
        for grp, lab in ((("H", "G", "I"), "helix"), (("E", "B"), "strand"),
                         (("T", "S", "C", "-", ""), "loop or coil")):
            vs = [v for _, _, v, s in fp if s in grp]
            if len(vs) >= 30:
                print("   %-14s n=%-5d mean cost %+.3f kcal/mol, median %+.3f"
                      % (lab, len(vs), st.mean(vs), st.median(vs)))

    # ------------------------------------------------------------------ 4
    print("\n" + "=" * 96)
    print("4  TRANSITION AGAINST TRANSVERSION, WITHIN A SOURCE RESIDUE\n")
    hs = parse_humsavar(os.path.expanduser(a.humsavar_out))
    if not hs:
        print("could not read the humsavar run at %s" % a.humsavar_out)
    else:
        pairs = [("D", "D>N", "D>H"), ("N", "N>D", "N>H"),
                 ("E", "E>K", "E>Q"), ("K", "K>E", "K>Q")]
        print("%-8s %-22s %-22s %s" % ("source", "transition", "transversion", "ratio"))
        ratios = []
        for src, ts, tv in pairs:
            if ts not in hs or tv not in hs:
                continue
            nts, ntv = sum(hs[ts]), sum(hs[tv])
            ratios.append(nts / ntv if ntv else float("nan"))
            print("%-8s %-10s n=%-9d %-10s n=%-9d %.2f" % (src, ts, nts, tv, ntv, ratios[-1]))
        if ratios:
            print("\nmedian transition to transversion ratio across the four source residues: %.2f"
                  % st.median(ratios))
            out_rows.append(("transition", "median_ratio", "%.2f" % st.median(ratios)))

    with open(os.path.expanduser(a.out), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["analysis", "key", "value"])
        for r in out_rows:
            w.writerow(r)
    print("\nwrote", os.path.expanduser(a.out))


if __name__ == "__main__":
    main()
