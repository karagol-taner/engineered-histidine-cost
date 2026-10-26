#!/usr/bin/env python3
"""Three checks on the pooled ordering.

Ranks the substitutions inside each protein separately, resamples proteins rather
than positions for the confidence intervals, and stratifies positions by a
leave-one-out constraint score.

    python3 ladder_robustness.py --help
"""
import argparse, csv, glob, gzip, importlib.util, os, sys
import numpy as np

NBIN = 2000
BOOT = 1000
TOL = 0.005


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


def cliff_hist(hx, hb):
    """Cliff's delta of a value histogram against a baseline histogram."""
    nx, nb = hx.sum(), hb.sum()
    if nx == 0 or nb == 0:
        return float("nan")
    cum = np.cumsum(hb)
    lt = np.concatenate(([0.0], cum[:-1]))      # baseline strictly below each bin
    gt = nb - cum                               # baseline strictly above each bin
    return float((hx * (gt - lt)).sum() / (nx * nb))


def read_exact(path):
    """The exact pooled Cliff deltas, for checking the binning."""
    out, inblock = {}, False
    if not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith("--- 12 core edges"):
            inblock = True
            continue
        if inblock and line.startswith("---"):
            break
        f = line.split()
        if inblock and len(f) > 6 and ">" in f[0]:
            try:
                out[f[0]] = float(f[6])
            except ValueError:
                pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="",
                    help="directory holding the method's own "
                         "modules; omit to use the bundled "
                         "codes/ladder_sets.py")
    ap.add_argument("--cohort", nargs="*",
                    default=["../../data/alphamissense_soluble.csv.gz"])
    ap.add_argument("--label", default="soluble")
    ap.add_argument("--exact", default="")
    ap.add_argument("--min-per-edge", type=int, default=5)
    ap.add_argument("--min-edges", type=int, default=12)
    ap.add_argument("--deciles", type=int, default=10)
    ap.add_argument("--boot", type=int, default=BOOT)
    ap.add_argument("--seed", type=int, default=20261007)
    ap.add_argument("--outdir", default="../../data")
    a = ap.parse_args()

    lec = load(tpath(a.targets, "ladder_edge_cost.py"), "lec")
    ca = load(tpath(a.targets, "codon_audit.py"), "ca")
    CORE, BASICDIR, SCHEMAS = lec.CORE, lec.BASICDIR, lec.SCHEMAS
    assert all(ca.routes(p[0], p[1]) for p in CORE), "a core edge is not single-base reachable"
    EDGE = [p[0] + ">" + p[1] for p in CORE]
    EIX = {e: i for i, e in enumerate(EDGE)}

    print("ladder_robustness.py")
    print("cohort = %s, %d bins, %d bootstrap resamples over proteins\n"
          % (a.label, NBIN, a.boot))

    paths = []
    for p in a.cohort:
        paths.extend(sorted(glob.glob(os.path.expanduser(p))) or [os.path.expanduser(p)])

    # ------------------------------------------------------------ one pass
    prot_ix, pos_ix = {}, {}
    psum, pcnt = [], []                      # per position, for the constraint score
    base_p, base_pos, base_v = [], [], []    # baseline rows
    edge_e, edge_p, edge_pos, edge_v = [], [], [], []
    total = 0
    for path in paths:
        op = gzip.open(path, "rt") if path.endswith(".gz") else open(path, newline="")
        with op as f:
            rd = csv.reader(f)
            header = [h.lstrip("# ").strip() for h in next(rd)]
            sc = next((c for c in SCHEMAS if c["aa1"] in header and c["score"] in header), None)
            if sc is None:
                sys.exit("unrecognised header in %s" % path)
            iacc = next((header.index(k) for k in ("acc", "Uniprot ACC") if k in header), None)
            ipos = header.index("position") if "position" in header else None
            if iacc is None or ipos is None:
                sys.exit("%s carries no accession/position columns" % path)
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
                total += 1
                acc = r[iacc]
                pi = prot_ix.get(acc)
                if pi is None:
                    pi = prot_ix[acc] = len(prot_ix)
                key = (pi, r[ipos])
                qi = pos_ix.get(key)
                if qi is None:
                    qi = pos_ix[key] = len(pos_ix)
                    psum.append(0.0); pcnt.append(0)
                psum[qi] += v; pcnt[qi] += 1
                snv = (r[isnv].strip().lower() == "y") if isnv is not None else True
                if snv:
                    base_p.append(pi); base_pos.append(qi); base_v.append(v)
                e = EIX.get(r[i1] + ">" + r[i2])
                if e is not None:
                    edge_e.append(e); edge_p.append(pi); edge_pos.append(qi); edge_v.append(v)
    if total == 0:
        sys.exit("no scored rows parsed")
    NP = len(prot_ix)
    psum = np.array(psum, dtype=np.float64); pcnt = np.array(pcnt, dtype=np.int32)
    base_p = np.array(base_p, dtype=np.int32); base_pos = np.array(base_pos, dtype=np.int32)
    base_v = np.array(base_v, dtype=np.float32)
    edge_e = np.array(edge_e, dtype=np.int8); edge_p = np.array(edge_p, dtype=np.int32)
    edge_pos = np.array(edge_pos, dtype=np.int32); edge_v = np.array(edge_v, dtype=np.float32)
    print("parsed %d rows, %d proteins, %d positions, %d baseline, %d on a core edge"
          % (total, NP, len(pcnt), len(base_v), len(edge_v)))
    print("substitutions scored per position: median %d, min %d, max %d"
          % (int(np.median(pcnt)), pcnt.min(), pcnt.max()))

    # ------------------------------------------------- per-protein histograms
    lo, hi = 0.0, 1.0
    b_of = lambda v: np.clip(((v - lo) / (hi - lo) * NBIN).astype(np.int32), 0, NBIN - 1)
    bb, eb = b_of(base_v), b_of(edge_v)
    HB = np.bincount(base_p.astype(np.int64) * NBIN + bb,
                     minlength=NP * NBIN).reshape(NP, NBIN).astype(np.float32)
    HE = np.bincount((edge_e.astype(np.int64) * NP + edge_p) * NBIN + eb,
                     minlength=12 * NP * NBIN).reshape(12, NP, NBIN).astype(np.float32)

    pooled_b = HB.sum(axis=0)
    pooled = {EDGE[i]: cliff_hist(HE[i].sum(axis=0), pooled_b) for i in range(12)}

    exact = read_exact(os.path.expanduser(a.exact)) if a.exact else {}
    if exact:
        worst = max(abs(pooled[e] - exact[e]) for e in EDGE if e in exact)
        print("\nbinning check against the exact pooled values: largest difference %.4f "
              "over %d edges" % (worst, len(exact)))
        if worst > TOL:
            sys.exit("FATAL: the %d-bin histogram disagrees with the exact Cliff delta by "
                     "%.4f, above the %.3f tolerance. Refusing to bootstrap on it."
                     % (NBIN, worst, TOL))
    else:
        print("\nno exact run supplied, so the binning is unchecked")

    order = sorted(EDGE, key=lambda e: -pooled[e])
    print("\npooled ordering (histogram): " + ", ".join(order))

    # ======================================================= 1 within-protein
    print("\n" + "=" * 100)
    print("1  WITHIN-PROTEIN REPLICATION")
    print("a protein enters with at least %d positions for each of at least %d "
          "substitutions\n" % (a.min_per_edge, a.min_edges))
    per_n = HE.sum(axis=2)                      # (12, NP)
    rhos, dh_last, nh_beats_dh, used = [], 0, 0, []
    perprot, percell = [], []
    for p in range(NP):
        ok = [i for i in range(12) if per_n[i, p] >= a.min_per_edge]
        if len(ok) < a.min_edges or HB[p].sum() < 200:
            continue
        loc = {EDGE[i]: cliff_hist(HE[i, p], HB[p]) for i in ok}
        rhos.append(spearman([loc[EDGE[i]] for i in ok], [pooled[EDGE[i]] for i in ok]))
        used.append(p)
        o = sorted(loc, key=lambda e: -loc[e])
        if o[-1] == "D>H":
            dh_last += 1
        if "N>H" in loc and "D>H" in loc and loc["N>H"] > loc["D>H"]:
            nh_beats_dh += 1
        perprot.append((p, int(HB[p].sum()), len(ok), rhos[-1],
                        o.index("D>H") + 1 if "D>H" in loc else "",
                        (loc["N>H"] - loc["D>H"]) if ("N>H" in loc and "D>H" in loc) else ""))
        for i in ok:
            percell.append((p, EDGE[i], int(per_n[i, p]), loc[EDGE[i]]))
    if len(used) < 10:
        print("only %d proteins qualify, which is too few to report" % len(used))
    else:
        r = np.array(rhos)
        print("proteins qualifying: %d of %d" % (len(used), NP))
        print("Spearman rho against the pooled ordering, per protein:")
        print("   median %+.3f, quartiles %+.3f to %+.3f, range %+.3f to %+.3f"
              % (np.median(r), np.percentile(r, 25), np.percentile(r, 75), r.min(), r.max()))
        print("   positive in %d of %d proteins (%.1f%%)"
              % ((r > 0).sum(), len(r), 100.0 * (r > 0).sum() / len(r)))
        print("   above +0.5 in %d (%.1f%%)" % ((r > 0.5).sum(), 100.0 * (r > 0.5).sum() / len(r)))
        print("D>H is the worst of the substitutions present in %d of %d proteins (%.1f%%)"
              % (dh_last, len(used), 100.0 * dh_last / len(used)))
        print("N>H is better tolerated than D>H in %d of %d proteins (%.1f%%)"
              % (nh_beats_dh, len(used), 100.0 * nh_beats_dh / len(used)))

    # ========================================================= 2 bootstrap
    print("\n" + "=" * 100)
    print("2  CLUSTER BOOTSTRAP OVER PROTEINS, %d resamples" % a.boot)
    rng = np.random.default_rng(a.seed)
    HEflat = np.ascontiguousarray(HE.transpose(0, 2, 1).reshape(12 * NBIN, NP))
    HBt = np.ascontiguousarray(HB.T)                       # (NBIN, NP)
    bd = np.zeros((a.boot, 12)); brank = np.zeros((a.boot, 12), dtype=np.int16)
    CH = 100                                               # resamples per block
    for s0 in range(0, a.boot, CH):
        nb_blk = min(CH, a.boot - s0)
        C = np.zeros((NP, nb_blk), dtype=np.float32)
        for j in range(nb_blk):
            C[:, j] = np.bincount(rng.integers(0, NP, NP), minlength=NP)
        hb = HBt @ C                                       # (NBIN, blk)
        cum = np.cumsum(hb, axis=0)
        nbv = cum[-1]
        w = (nbv[None, :] - cum) - np.vstack([np.zeros((1, nb_blk), np.float32), cum[:-1]])
        he = (HEflat @ C).reshape(12, NBIN, nb_blk)        # (12, NBIN, blk)
        nx = he.sum(axis=1)                                # (12, blk)
        num = np.einsum('enb,nb->eb', he, w)
        d = np.where(nx > 0, num / np.maximum(nx * nbv[None, :], 1.0), np.nan)
        bd[s0:s0 + nb_blk] = d.T
        brank[s0:s0 + nb_blk] = ((-d).argsort(axis=0).argsort(axis=0) + 1).T
    print("%-7s %-9s %-20s %-8s %s" % ("edge", "pooled", "95% interval", "rank", "rank 95%"))
    for e in order:
        i = EIX[e]
        v = bd[:, i]
        pr = order.index(e) + 1
        print("%-7s %+-9.3f %+.3f to %+.3f      %-8d %d to %d"
              % (e, pooled[e], np.percentile(v, 2.5), np.percentile(v, 97.5), pr,
                 int(np.percentile(brank[:, i], 2.5)), int(np.percentile(brank[:, i], 97.5))))
    dh, nh = EIX["D>H"], EIX["N>H"]
    print("\nD>H is last of twelve in %.1f%% of resamples" % (100.0 * (brank[:, dh] == 12).mean()))
    print("D>H is in the bottom two in %.1f%% of resamples" % (100.0 * (brank[:, dh] >= 11).mean()))
    print("N>H is better tolerated than D>H in %.1f%% of resamples"
          % (100.0 * (bd[:, nh] > bd[:, dh]).mean()))
    print("the N>H minus D>H difference is %+.3f, 95%% interval %+.3f to %+.3f"
          % (pooled["N>H"] - pooled["D>H"],
             np.percentile(bd[:, nh] - bd[:, dh], 2.5),
             np.percentile(bd[:, nh] - bd[:, dh], 97.5)))

    # ============================================== 3 constraint stratification
    print("\n" + "=" * 100)
    print("3  LEAVE-ONE-OUT POSITIONAL CONSTRAINT")
    print("deciles are cut on the baseline rows\n")
    keep = pcnt > 1
    if keep.sum() < 1000:
        sys.exit("FATAL: too few positions carry more than one scored substitution")
    loo = lambda pos, v: (psum[pos] - v) / np.maximum(pcnt[pos] - 1, 1)
    cb = loo(base_pos, base_v.astype(np.float64))
    ce = loo(edge_pos, edge_v.astype(np.float64))
    cuts = np.percentile(cb, np.linspace(0, 100, a.deciles + 1))
    cuts[0] -= 1e-9; cuts[-1] += 1e-9
    db = np.clip(np.searchsorted(cuts, cb, side="right") - 1, 0, a.deciles - 1)
    de = np.clip(np.searchsorted(cuts, ce, side="right") - 1, 0, a.deciles - 1)
    print("constraint of a position, mean of the other substitutions: median %.3f, "
          "decile cuts %s" % (np.median(cb), " ".join("%.2f" % c for c in cuts[1:-1])))
    bbv, ebv = b_of(base_v), b_of(edge_v)
    print("\n%-7s %s" % ("edge", "  ".join("D%-2d" % (k + 1) for k in range(a.deciles))))
    strat = {}
    for e in order:
        i = EIX[e]
        row, sel_e = [], (edge_e == i)
        for k in range(a.deciles):
            m = sel_e & (de == k)
            hb = np.zeros(NBIN); np.add.at(hb, bbv[db == k], 1.0)
            hx = np.zeros(NBIN); np.add.at(hx, ebv[m], 1.0)
            row.append(cliff_hist(hx, hb) if m.sum() >= 30 else float("nan"))
        strat[e] = row
        print("%-7s %s" % (e, "  ".join(("%+.2f" % v) if v == v else "  . " for v in row)))
    print("\nrank of each substitution inside each decile (1 = best tolerated):")
    print("%-7s %s" % ("edge", "  ".join("D%-2d" % (k + 1) for k in range(a.deciles))))
    ranks = {e: [] for e in EDGE}
    for k in range(a.deciles):
        col = [(e, strat[e][k]) for e in EDGE if strat[e][k] == strat[e][k]]
        col.sort(key=lambda t: -t[1])
        pos = {e: j + 1 for j, (e, _) in enumerate(col)}
        for e in EDGE:
            ranks[e].append(pos.get(e, float("nan")))
    for e in order:
        print("%-7s %s" % (e, "  ".join(("%-3d" % r) if r == r else " .  " for r in ranks[e])))
    dh_ranks = [r for r in ranks["D>H"] if r == r]
    print("\nD>H is last in %d of the %d reportable deciles" % (sum(1 for r in dh_ranks if r == 12),
                                                               len(dh_ranks)))
    print("D>H's Cliff delta runs from %+.3f in the least constrained decile to %+.3f in the most"
          % (strat["D>H"][0], strat["D>H"][-1]))
    print("N>H runs from %+.3f to %+.3f over the same range" % (strat["N>H"][0], strat["N>H"][-1]))
    rr = [spearman([strat[e][k] for e in EDGE if strat[e][k] == strat[e][k]],
                   [pooled[e] for e in EDGE if strat[e][k] == strat[e][k]])
          for k in range(a.deciles)]
    print("Spearman rho of each decile's ordering against the pooled ordering: %s"
          % ", ".join("%+.2f" % v for v in rr))

    pc = os.path.join(os.path.expanduser(a.outdir), "ladder_percell_%s.csv" % a.label)
    with open(pc, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["protein_index", "edge", "n", "cliff_within_protein"])
        for row in percell:
            w.writerow([row[0], row[1], row[2], "%.4f" % row[3]])
    print("wrote", pc)

    wp = os.path.join(os.path.expanduser(a.outdir), "ladder_within_protein_%s.csv" % a.label)
    with open(wp, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["protein_index", "baseline_n", "edges_scored", "spearman_rho_vs_pooled",
                    "rank_of_DtoH", "NtoH_minus_DtoH"])
        for row in perprot:
            w.writerow(list(row[:3]) + ["%.4f" % row[3], row[4],
                                        ("%.4f" % row[5]) if row[5] != "" else ""])
    print("\nwrote", wp)

    out = os.path.join(os.path.expanduser(a.outdir), "ladder_robustness_%s.csv" % a.label)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["edge", "direction", "cohort", "pooled_cliff", "boot_lo", "boot_hi",
                    "pooled_rank", "rank_lo", "rank_hi"]
                   + ["constraint_decile_%d" % (k + 1) for k in range(a.deciles)])
        for p in CORE:
            e = p[0] + ">" + p[1]
            i = EIX[e]
            w.writerow([e, "basic" if p in BASICDIR else "acidic", a.label,
                        "%.3f" % pooled[e],
                        "%.3f" % np.percentile(bd[:, i], 2.5),
                        "%.3f" % np.percentile(bd[:, i], 97.5),
                        order.index(e) + 1,
                        int(np.percentile(brank[:, i], 2.5)),
                        int(np.percentile(brank[:, i], 97.5))]
                       + [("%.3f" % v) if v == v else "" for v in strat[e]])
    print("\nwrote", out)


if __name__ == "__main__":
    main()
