#!/usr/bin/env python3
import argparse, collections, csv, os, sys

import importlib.util


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--megascale", default="../../data/megascale")
    ap.add_argument("--core", default="../../data/ladder_core_fullprec.csv")
    ap.add_argument("--targets", default="")
    ap.add_argument("--out", default="../../data/ladder_source_null.csv")
    ap.add_argument("--min-n", type=int, default=20)
    a = ap.parse_args()

    lec = _sets(a.targets, "ladder_edge_cost.py")
    ca = _sets(a.targets, "codon_audit.py")
    CORE = [tuple(p) for p in lec.CORE]
    snv_ok = lambda p: bool(ca.routes(p[0], p[1]))

    raw = LM.read_megascale(a.megascale)
    flip = (sum(1 for r in raw if r[4] < 0) / len(raw)) > 0.5
    rows = [(d, p, a1, a2, (-v if flip else v)) for d, p, a1, a2, v, _ in raw]
    print("ladder_source_null.py: the source-residue null against the twelve")
    print("cDNA-proteolysis single substitutions %d, costs mapped with higher = "
          "less stable" % len(rows))

    base = [v for _, _, a1, a2, v in rows if snv_ok((a1, a2))]
    cliff, nb = LM.cliff_factory(base)
    per = collections.defaultdict(list)
    for _, _, a1, a2, v in rows:
        per[(a1, a2)].append(v)
    cost = {p: (len(xs), cliff(xs)) for p, xs in per.items() if len(xs) >= a.min_n}
    print("baseline every single-nucleotide-reachable substitution, n=%d; "
          "%d substitution classes with at least %d measurements"
          % (nb, len(cost), a.min_n))


    pub = {}
    for r in csv.DictReader(open(a.core, encoding="utf-8")):
        s, d = r["edge"].split(">")
        pub[(s, d)] = float(r["cliff_megascale"])
    worst = max(abs(cost[p][1] - pub[p]) for p in CORE)
    print("gate: the twelve against %s, largest difference %.6f"
          % (os.path.basename(a.core), worst))
    if worst > 1e-3:
        sys.exit("FATAL: the recomputed ordering does not match %s" % a.core)

    def source_effect(s, exclude):
        sub = [(n, c) for (a1, a2), (n, c) in cost.items() if a1 == s and a2 != exclude]
        N = sum(n for n, _ in sub)
        return sum(n * c for n, c in sub) / N, N, len(sub)

    out = []
    for s, d in CORE:
        n, c = cost[(s, d)]
        eff, N, k = source_effect(s, d)
        out.append({"pair": "%s>%s" % (s, d), "n": n, "cliff": c,
                    "source_effect": eff, "source_n": N, "source_classes": k,
                    "residual": c - eff})

    obs = [r["cliff"] for r in out]
    src = [r["source_effect"] for r in out]

    try:
        from Bio.Align import substitution_matrices
        B = substitution_matrices.load("BLOSUM62")
        bl = [float(B[r["pair"][0], r["pair"][2]]) for r in out]
        for r, v in zip(out, bl):
            r["blosum62"] = v
        have_bl = True
    except Exception as e:
        print("BLOSUM62 unavailable (%s), that null is not reported" % e)
        bl, have_bl = None, False

    print("\n%-6s %7s %9s %9s %9s" % ("pair", "n", "cliff", "source", "residual"))
    for r in sorted(out, key=lambda r: -r["cliff"]):
        print("%-6s %7d %+9.3f %+9.3f %+9.3f"
              % (r["pair"], r["n"], r["cliff"], r["source_effect"], r["residual"]))

    print("\nrank agreement with the observed ordering of the twelve")
    print("   source-residue only        rho %+.3f" % LM.spearman(obs, src))
    if have_bl:
        print("   BLOSUM62                   rho %+.3f" % LM.spearman(obs, bl))

    print("\nN>H against D>H")
    nh = next(x for x in out if x["pair"] == "N>H")
    dh = next(x for x in out if x["pair"] == "D>H")
    for r in (nh, dh):
        print("   %s   cliff %+.3f   source %s>* %+.3f   residual %+.3f"
              % (r["pair"], r["cliff"], r["pair"][0], r["source_effect"], r["residual"]))
    gap, sgap = nh["cliff"] - dh["cliff"], nh["source_effect"] - dh["source_effect"]
    print("   observed gap %+.3f, of which the source residues explain %+.3f (%.0f%%),"
          % (gap, sgap, 100.0 * sgap / gap))
    print("   leaving %+.3f once the source residue is removed" % (gap - sgap))

    print("\ndirectional gap between reciprocal pairs")
    seen, asym = set(), []
    for s1, d1 in CORE:
        if (s1, d1) in seen or (d1, s1) in seen or (d1, s1) not in cost:
            continue
        seen.add((s1, d1))
        f, b = cost[(s1, d1)][1], cost[(d1, s1)][1]
        asym.append((s1, d1, f, b, f - b))
    for s1, d1, f, b, g in sorted(asym, key=lambda t: -abs(t[4])):
        print("   %s>%s %+.3f   %s>%s %+.3f   direction %+.3f" % (s1, d1, f, d1, s1, b, g))
    spread = max(cost[p][1] for p in CORE) - min(cost[p][1] for p in CORE)
    mean_abs = sum(abs(t[4]) for t in asym) / len(asym)
    print("   %d reciprocal pairs, mean directional gap %.3f, which is %.0f%% of "
          "the %.3f spread across the twelve"
          % (len(asym), mean_abs, 100.0 * mean_abs / spread, spread))
    print("   BLOSUM62 and Grantham are symmetric, so each of these gaps is zero "
          "under either matrix")

    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        cols = ["pair", "n", "cliff", "source_effect", "source_n", "source_classes",
                "residual"] + (["blosum62"] if have_bl else [])
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in sorted(out, key=lambda r: -r["cliff"]):
            w.writerow({k: r[k] for k in cols})
    print("\nwritten: %s" % a.out)


if __name__ == "__main__":
    main()
