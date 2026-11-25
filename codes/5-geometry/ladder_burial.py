#!/usr/bin/env python3
import argparse, csv, importlib.util, os, sys


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


def identifiers(path):
    with open(path, encoding="utf-8") as fh:
        head = fh.readline()
    if "pdb" in [c.strip().lower() for c in head.split(",")]:
        with open(path, newline="", encoding="utf-8") as fh:
            return [r["pdb"].strip().upper() for r in csv.DictReader(fh)
                    if r.get("pdb", "").strip()]
    return [l.strip().upper() for l in open(path, encoding="utf-8") if l.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--megascale", default="../../data/megascale")
    ap.add_argument("--asa", default="../../data/ladder_burial_asa.csv")
    ap.add_argument("--ids", default="../../data/ladder_geometry_structures.csv")
    ap.add_argument("--core", default="../../data/ladder_core_fullprec.csv")
    ap.add_argument("--targets", default="")
    ap.add_argument("--buried-asa", type=float, default=25.0)
    ap.add_argument("--min-identity", type=float, default=0.8)
    ap.add_argument("--out", default="../../data/ladder_burial.csv")
    a = ap.parse_args()

    lec = _sets(a.targets, "ladder_edge_cost.py")
    ca = _sets(a.targets, "codon_audit.py")
    CORE = [tuple(p) for p in lec.CORE]
    snv_ok = lambda p: bool(ca.routes(p[0], p[1]))

    raw = LM.read_megascale(a.megascale)
    frac_neg = sum(1 for r in raw if r[4] < 0) / len(raw)
    flip = frac_neg > 0.5
    print("ladder_burial.py: does burial produce the ordering of the twelve?")
    print("megascale single substitutions %d, %.1f%% negative, so %s is "
          "destabilising, mapped to cost with higher = less stable"
          % (len(raw), 100 * frac_neg, "negative" if flip else "positive"))
    cost = [(d, p, a1, a2, (-v if flip else v), s) for d, p, a1, a2, v, s in raw]


    pooled, _ = LM.rank_edges([(d, p, a1, a2, v) for d, p, a1, a2, v, _ in cost], snv_ok, CORE)
    pub = {}
    for r in csv.DictReader(open(a.core, encoding="utf-8")):
        s_, d_ = r["edge"].split(">")
        pub[(s_, d_)] = float(r["cliff_megascale"])
    worst = max(abs(pooled[p][0] - pub[p]) for p in CORE)
    print("   gate: the twelve against %s, largest difference %.6f"
          % (os.path.basename(a.core), worst))
    if worst > 1e-3:
        sys.exit("FATAL: the recomputed ordering does not match %s" % a.core)


    seqs = LM.chain_sequences(a.asa)
    want = set(identifiers(a.ids))
    dom_seq = {}
    for d, p, a1, a2, v, s in cost:
        if d in want and s and d not in dom_seq:
            dom_seq[d] = s
    placed, skipped = {}, []
    for d, s in dom_seq.items():
        cands = [(k, vv) for k, vv in seqs.items() if k[0] == d]
        best = (None, 0.0, None)
        for k, (cs, asas) in cands:
            o, idn = LM.best_offset(s, cs)
            if o is not None and idn > best[1]:
                best = (k, idn, o)
        if best[0] and best[1] >= a.min_identity:
            placed[d] = (best[0], best[2], best[1])
        else:
            skipped.append((d, round(best[1], 2)))
    print("   natural domains %d, placed on a chain at >= %.0f%% identity %d, "
          "unplaced %d" % (len(dom_seq), 100 * a.min_identity, len(placed), len(skipped)))

    placed_rows, wrong = [], 0
    for d, p, a1, a2, v, _s in cost:
        if d not in placed:
            continue
        key, off, _ = placed[d]
        cs, asas = seqs[key]
        i = off + p - 1
        if i < 0 or i >= len(asas):
            continue
        if cs[i] != a1:
            wrong += 1
            continue
        placed_rows.append((d, p, a1, a2, v, asas[i]))
    print("   measurements placed on a residue %d, dropped for a wild-type "
          "mismatch %d" % (len(placed_rows), wrong))


    areas = sorted(r[5] for r in placed_rows)

    def pct(q):
        return areas[min(len(areas) - 1, int(q * len(areas)))]

    cuts = [("%.0f A^2 (ladder_ddg.py)" % a.buried_asa, a.buried_asa)]
    for q in (0.20, 0.30, 0.40, 0.50, 0.60):
        cuts.append(("%.0f%% most buried, %.1f A^2" % (100 * q, pct(q)), pct(q)))

    print("\nthe ordering under each burial threshold, with the baseline of the "
          "stratum under it")
    print("   %-28s %7s %8s %7s %9s" % ("threshold", "buried", "exposed",
                                        "edges", "rho b vs e"))
    rows_out, best = [], None
    for label, thr in cuts:
        bur = [(d, p, x, y, v) for d, p, x, y, v, s_ in placed_rows if s_ < thr]
        exp = [(d, p, x, y, v) for d, p, x, y, v, s_ in placed_rows if s_ >= thr]
        rb, _ = LM.rank_edges(bur, snv_ok, CORE)
        re_, _ = LM.rank_edges(exp, snv_ok, CORE)
        if not rb or not re_:
            print("   %-28s %7d %8d   too few to rank" % (label, len(bur), len(exp)))
            continue
        common = [p for p in CORE if p in rb and p in re_]
        rho_be = LM.spearman([rb[p][0] for p in common], [re_[p][0] for p in common])
        print("   %-28s %7d %8d %7d %+9.3f"
              % (label, len(bur), len(exp), len(common), rho_be))
        rows_out.append((label, thr, len(bur), len(exp), len(common), rho_be))
        if best is None or len(common) > best[4]:
            best = (label, thr, len(bur), len(exp), len(common), rho_be, rb, re_, common)

    label, thr, nb_, ne_, nc, rho_be, rb, re_, common = best
    print("\n--- the threshold ranking the most edges in both strata: %s ---" % label)
    for name, r in (("buried", rb), ("exposed", re_)):
        print("   %s" % name)
        for p in sorted(r, key=lambda p: -r[p][0]):
            print("      %s>%s  n=%5d  cliff %+.3f" % (p[0], p[1], r[p][1], r[p][0]))
    pl = [pooled[p][0] for p in common]
    print("\n   buried against exposed   rho %+.3f over %d edges" % (rho_be, nc))
    print("   buried against pooled    rho %+.3f"
          % LM.spearman([rb[p][0] for p in common], pl))
    print("   exposed against pooled   rho %+.3f"
          % LM.spearman([re_[p][0] for p in common], pl))

    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["threshold", "area", "n_buried", "n_exposed", "edges", "rho_buried_vs_exposed"])
        for row in rows_out:
            w.writerow([row[0], round(row[1], 2), row[2], row[3], row[4], round(row[5], 4)])
        w.writerow([])
        w.writerow(["pair", "n_buried", "cliff_buried", "n_exposed", "cliff_exposed", "cliff_pooled"])
        for p in common:
            w.writerow(["%s>%s" % p, rb[p][1], round(rb[p][0], 4),
                        re_[p][1], round(re_[p][0], 4), round(pooled[p][0], 4)])
    print("\nwritten: %s" % a.out)


if __name__ == "__main__":
    main()
