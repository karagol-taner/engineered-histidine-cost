#!/usr/bin/env python3
"""Functional stratification of the substitution set.

Splits a variant cohort into positions UniProt annotates as an active site, a
binding site or a site and the rest, and scores each substitution against the
baseline computed inside its own stratum.

    python3 ladder_functional.py --help
"""
import argparse, bisect, csv, glob, gzip, importlib.util, json, os
import statistics as st, sys, time, urllib.request

FEATURE_TYPES = {"Active site": "ACT_SITE", "Binding site": "BINDING", "Site": "SITE"}
BATCH = 90


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
    base = sorted(base)
    nb = len(base)

    def cliff(xs):
        gt = lt = 0
        for x in xs:
            lt += bisect.bisect_left(base, x)
            gt += nb - bisect.bisect_right(base, x)
        return (gt - lt) / (len(xs) * nb)
    return cliff, nb


def wilson(k, n):
    if n == 0:
        return (float('nan'), float('nan'))
    z, p = 1.959964, k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5)
    return ((c - h) / d * 100, (c + h) / d * 100)


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


# --------------------------------------------------------------- the cohort
def read_cohort(paths, wanted, SCHEMAS):
    """One pass over any AlphaMissense cohort in either schema."""
    base_all, rows, accs = [], {p: [] for p in wanted}, set()
    acc_keys = ('acc', 'Uniprot ACC')
    for path in paths:
        op = gzip.open(path, 'rt') if path.endswith('.gz') else open(path, newline='')
        with op as f:
            rd = csv.reader(f)
            header = [h.lstrip('# ').strip() for h in next(rd)]
            sc = None
            for cand in SCHEMAS:
                if cand['aa1'] in header and cand['score'] in header:
                    sc = cand
                    break
            if sc is None:
                sys.exit('unrecognised header in %s: %s' % (path, header[:8]))
            iacc = next((header.index(k) for k in acc_keys if k in header), None)
            ipos = header.index('position') if 'position' in header else None
            if iacc is None or ipos is None:
                sys.exit('%s carries no accession/position columns, so positions '
                         'cannot be annotated' % path)
            i1, i2 = header.index(sc['aa1']), header.index(sc['aa2'])
            isc = header.index(sc['score'])
            isnv = header.index(sc['snv']) if sc['snv'] in header else None
            ipc = header.index(sc['pclass']) if sc['pclass'] in header else None
            for r in rd:
                if len(r) <= isc:
                    continue
                try:
                    v = float(r[isc])
                except ValueError:
                    continue
                snv = (r[isnv].strip().lower() == 'y') if isnv is not None else True
                acc, pos = r[iacc], r[ipos]
                if snv:
                    base_all.append((acc, pos, v))
                    accs.add(acc)
                k = (r[i1], r[i2])
                if k in rows:
                    rows[k].append((acc, pos, v, r[ipc] if ipc is not None else '', r[i1]))
                    accs.add(acc)
    if not base_all:
        sys.exit('no single-nucleotide-reachable rows parsed; check the is_snv column')
    return base_all, rows, accs


# ------------------------------------------------------------ uniprot features
def fetch_features(accs, cache):
    """ACT_SITE, BINDING and SITE positions per accession, plus the sequence for
    the position check."""
    have = {}
    if os.path.exists(cache):
        with open(cache) as fh:
            have = json.load(fh)
    todo = [a for a in sorted(accs) if a not in have]
    if todo:
        print("fetching UniProt features for %d accessions (%d cached), %d per request"
              % (len(todo), len(have), BATCH))
    for s in range(0, len(todo), BATCH):
        chunk = todo[s:s + BATCH]
        url = ("https://rest.uniprot.org/uniprotkb/accessions?accessions=" + ",".join(chunk)
               + "&fields=accession,ft_act_site,ft_binding,ft_site,sequence&format=json")
        for attempt in (1, 2, 3):
            try:
                with urllib.request.urlopen(url, timeout=120) as r:
                    d = json.load(r)
                break
            except Exception as e:
                if attempt == 3:
                    sys.exit("FATAL: UniProt fetch failed after 3 attempts (%s). A partial "
                             "annotation would silently mislabel positions as unannotated, "
                             "so the run stops." % e)
                time.sleep(3 * attempt)
        got = set()
        for e in d.get("results", []):
            acc = e.get("primaryAccession")
            got.add(acc)
            feats = {}
            for ft in e.get("features", []):
                t = FEATURE_TYPES.get(ft.get("type"))
                if not t:
                    continue
                loc = ft.get("location", {})
                a = (loc.get("start") or {}).get("value")
                b = (loc.get("end") or {}).get("value")
                if a is None or b is None:
                    continue
                for p in range(int(a), int(b) + 1):
                    feats.setdefault(str(p), []).append(t)
            have[acc] = {"feat": feats, "seq": (e.get("sequence") or {}).get("value", "")}
        for acc in chunk:                      # obsolete or demerged accessions
            if acc not in got:
                have[acc] = {"feat": {}, "seq": "", "missing": True}
        with open(cache, "w") as fh:
            json.dump(have, fh)
        print("   %d/%d" % (min(s + BATCH, len(todo)), len(todo)))
        time.sleep(0.3)
    miss = [a for a in accs if have.get(a, {}).get("missing")]
    if miss:
        print("   %d accessions returned nothing from UniProt and are treated as "
              "unannotated: %s" % (len(miss), ', '.join(miss[:8]) + (' ...' if len(miss) > 8 else '')))
    return have


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="",
                    help="directory holding the method's own "
                         "modules; omit to use the bundled "
                         "codes/ladder_sets.py")
    ap.add_argument("--cohort", nargs="*",
                    default=["../../data/alphamissense_soluble.csv.gz"])
    ap.add_argument("--label", default="soluble")
    ap.add_argument("--fireprot", default="../../data/fireprot")
    ap.add_argument("--cache", default="../../data/uniprot_functional_features.json")
    ap.add_argument("--window", type=int, default=0,
                    help="treat a position within N residues of a feature as annotated")
    ap.add_argument("--minn", type=int, default=20)
    ap.add_argument("--out", default="")
    ap.add_argument("--skip-fireprot", action="store_true")
    a = ap.parse_args()
    out = a.out or "../../data/ladder_functional_%s_w%d.csv" % (a.label, a.window)

    lec = load(tpath(a.targets, "ladder_edge_cost.py"), "lec")
    ca = load(tpath(a.targets, "codon_audit.py"), "ca")
    CORE, ARG, EXT, BASICDIR = lec.CORE, lec.ARG, lec.EXT, lec.BASICDIR
    wanted = set(CORE) | set(ARG) | set(EXT)
    assert all(ca.routes(p[0], p[1]) for p in CORE), "a core edge is not single-base reachable"

    print("ladder_functional.py")
    print("cohort = %s, annotation window = +/-%d residue(s)" % (a.label, a.window))

    paths = []
    for p in a.cohort:
        paths.extend(sorted(glob.glob(os.path.expanduser(p))) or [os.path.expanduser(p)])
    base_all, rows, accs = read_cohort(paths, wanted, lec.SCHEMAS)
    print("cohort: %d single-nucleotide baseline rows over %d accessions, %d file(s)"
          % (len(base_all), len(accs), len(paths)))

    feats = fetch_features(accs, os.path.expanduser(a.cache))

    checked = mismatch = 0
    for p, rs in rows.items():
        for acc, pos, _, _, a1 in rs:
            seq = feats.get(acc, {}).get("seq") or ""
            try:
                i = int(pos)
            except ValueError:
                continue
            if 1 <= i <= len(seq):
                checked += 1
                if seq[i - 1] != a1:
                    mismatch += 1
    rate = 100.0 * mismatch / checked if checked else float('nan')
    print("position check: %d of %d edge rows disagree with the UniProt residue (%.2f%%)"
          % (mismatch, checked, rate))
    if checked and rate > 2.0:
        sys.exit("FATAL: the cohort's numbering disagrees with UniProt at %.2f%% of "
                 "positions, so the annotation would land on the wrong residues." % rate)

    def annotated(acc, pos):
        f = feats.get(acc, {}).get("feat") or {}
        if not f:
            return False
        if not a.window:
            return pos in f
        try:
            i = int(pos)
        except ValueError:
            return False
        return any(str(j) in f for j in range(i - a.window, i + a.window + 1))

    base_ann = [v for acc, pos, v in base_all if annotated(acc, pos)]
    base_un = [v for acc, pos, v in base_all if not annotated(acc, pos)]
    nfeat = sum(len(v.get("feat") or {}) for v in feats.values())
    nprot = sum(1 for v in feats.values() if v.get("feat"))
    allv = [v for _, _, v in base_all]
    print("annotated residues in UniProt across the cohort: %d, in %d of %d proteins"
          % (nfeat, nprot, len(accs)))
    cliff_g, nb_g = cliff_factory(allv)
    print("baseline, whole cohort    n=%-9d mean=%.4f median=%.4f"
          % (nb_g, st.mean(allv), st.median(allv)))
    if len(base_ann) < 200 or len(base_un) < 200:
        sys.exit("FATAL: a stratum baseline is too small to compare against "
                 "(annotated n=%d, unannotated n=%d)." % (len(base_ann), len(base_un)))
    cliff_a, nb_a = cliff_factory(base_ann)
    cliff_u, nb_u = cliff_factory(base_un)
    print("baseline, annotated       n=%-9d mean=%.4f median=%.4f"
          % (nb_a, st.mean(base_ann), st.median(base_ann)))
    print("baseline, unannotated     n=%-9d mean=%.4f median=%.4f"
          % (nb_u, st.mean(base_un), st.median(base_un)))
    print("annotated share of the baseline: %.2f%%\n" % (100.0 * nb_a / nb_g))

    res = {}
    for p in wanted:
        rs = rows.get(p) or []
        ann = [(v, pc) for acc, pos, v, pc, _ in rs if annotated(acc, pos)]
        un = [(v, pc) for acc, pos, v, pc, _ in rs if not annotated(acc, pos)]
        rec = dict(n_ann=len(ann), n_un=len(un))
        if ann:
            va = [v for v, _ in ann]
            rec.update(mean_ann=st.mean(va), g_ann=cliff_g(va), m_ann=cliff_a(va),
                       path_ann=sum(1 for _, pc in ann if pc == 'likely_pathogenic'))
        if un:
            vu = [v for v, _ in un]
            rec.update(mean_un=st.mean(vu), g_un=cliff_g(vu), m_un=cliff_u(vu),
                       path_un=sum(1 for _, pc in un if pc == 'likely_pathogenic'))
        res[p] = rec

    def table(pairs, tag):
        print("=" * 112)
        print("--- %s ---" % tag)
        print("%-7s %-6s | %-8s %-8s %-8s %-8s | %-8s %-8s %-8s %-8s"
              % ("edge", "dir", "n_ann", "mean", "global", "matched",
                 "n_unann", "mean", "global", "matched"))
        got = []
        for p in pairs:
            r = res.get(p) or {}
            if r.get('n_ann', 0) < a.minn or r.get('n_un', 0) < a.minn:
                continue
            got.append(p)
            print("%-7s %-6s | %-8d %+-8.4f %+-8.3f %+-8.3f | %-8d %+-8.4f %+-8.3f %+-8.3f"
                  % (p[0] + '>' + p[1],
                     'basic' if p in BASICDIR else 'acidic' if (p[1], p[0]) in BASICDIR else '-',
                     r['n_ann'], r['mean_ann'], r['g_ann'], r['m_ann'],
                     r['n_un'], r['mean_un'], r['g_un'], r['m_un']))
        skipped = [p for p in pairs if p not in got]
        if skipped:
            print("below n=%d in a stratum and not reported: %s"
                  % (a.minn, ', '.join(x + '>' + y for x, y in skipped)))
        print()
        return got

    got = table(CORE, "12 core edges, by functional annotation")
    table(ARG, "5 arginine edges")
    table(EXT, "proposed extensions")

    if len(got) >= 3:
        print("=" * 112)
        print("rank within each stratum, matched baseline, 1 = best tolerated\n")
        ra = sorted(got, key=lambda p: -res[p]['m_ann'])
        ru = sorted(got, key=lambda p: -res[p]['m_un'])
        print("%-7s %-24s %-24s %s" % ("edge", "rank at annotated", "rank at unannotated", "shift"))
        for p in ra:
            print("%-7s %-24s %-24s %+d"
                  % (p[0] + '>' + p[1],
                     "%2d of %d  (%+.3f)" % (ra.index(p) + 1, len(ra), res[p]['m_ann']),
                     "%2d of %d  (%+.3f)" % (ru.index(p) + 1, len(ru), res[p]['m_un']),
                     ra.index(p) - ru.index(p)))
        rho = spearman([res[p]['m_ann'] for p in got], [res[p]['m_un'] for p in got])
        print("\nSpearman rho between the two strata's matched orderings: %+.3f over %d edges"
              % (rho, len(got)))

        print("\npathogenic fraction per stratum (AlphaMissense likely_pathogenic), "
              "Wilson 95%:")
        print("%-7s %-30s %-30s" % ("edge", "annotated", "unannotated"))
        for p in ra:
            r = res[p]
            la, ha = wilson(r['path_ann'], r['n_ann'])
            lu, hu = wilson(r['path_un'], r['n_un'])
            print("%-7s %-30s %-30s"
                  % (p[0] + '>' + p[1],
                     "%5.1f%%  (%.1f-%.1f, n=%d)" % (100.0 * r['path_ann'] / r['n_ann'], la, ha, r['n_ann']),
                     "%5.1f%%  (%.1f-%.1f, n=%d)" % (100.0 * r['path_un'] / r['n_un'], lu, hu, r['n_un'])))

    if not a.skip_fireprot:
        print("\n" + "=" * 112)
        print("MEASURED FOLDING COST BY FIREPROT'S OWN FUNCTIONAL FLAGS")
        print("flags: is_in_catalytic_pocket, is_essential\n")
        seen, fp = set(), []
        for path in sorted(glob.glob(os.path.join(a.fireprot, "*.csv"))):
            with open(path, newline="") as fh:
                for r in csv.DictReader(fh):
                    a1, a2 = r.get("wild_type"), r.get("mutation")
                    try:
                        v = float(r.get("ddG"))
                    except (TypeError, ValueError):
                        continue
                    if not a1 or not a2 or a1 == a2:
                        continue
                    key = (r.get("uniprot_id"), r.get("position"), a1, a2, round(v, 3))
                    if key in seen:
                        continue
                    seen.add(key)
                    flag = (str(r.get("is_in_catalytic_pocket")).strip().lower() in ('true', '1', 'yes')
                            or str(r.get("is_essential")).strip().lower() in ('true', '1', 'yes'))
                    fp.append((a1, a2, v, flag))
        if not fp:
            print("no FireProt rows read from %s" % a.fireprot)
        else:
            neg = sum(1 for _, _, v, _ in fp if v < 0) / len(fp)
            flip = neg > 0.5
            fp = [(x, y, (-v if flip else v), f) for x, y, v, f in fp]
            print("sign: %.1f%% of FireProt ddG values negative, so %s is destabilising"
                  % (100 * neg, "negative" if flip else "positive"))
            pro = [v for _, y, v, _ in fp if y == 'P']
            med = st.median([v for _, _, v, _ in fp])
            if len(pro) >= 30:
                if st.median(pro) <= med:
                    sys.exit("FATAL FireProt: X>P is not costlier than the median, so the "
                             "sign mapping is inverted.")
                print("sign check, proline: median %+.3f overall against %+.3f for X>P (n=%d)"
                      % (med, st.median(pro), len(pro)))
            nf = sum(1 for _, _, _, f in fp if f)
            print("%d of %d rows flagged catalytic-pocket or essential (%.1f%%)\n"
                  % (nf, len(fp), 100.0 * nf / len(fp)))
            print("%-7s | %-28s | %-28s" % ("edge", "flagged functional", "not flagged"))
            for p in CORE + [('Y', 'H'), ('N', 'K'), ('H', 'Q')]:
                fl = [v for x, y, v, f in fp if (x, y) == p and f]
                nn = [v for x, y, v, f in fp if (x, y) == p and not f]
                s = lambda xs: ("%+.3f kcal/mol (n=%d)" % (st.mean(xs), len(xs))) if xs else "no rows"
                print("%-7s | %-28s | %-28s" % (p[0] + '>' + p[1], s(fl), s(nn)))
            allf = [v for _, _, v, f in fp if f]
            alln = [v for _, _, v, f in fp if not f]
            print("\nALL substitutions, flagged functional   %+.3f kcal/mol (n=%d)"
                  % (st.mean(allf), len(allf)))
            print("ALL substitutions, not flagged          %+.3f kcal/mol (n=%d)"
                  % (st.mean(alln), len(alln)))
            print("This aggregate is the one well-powered folding comparison available:")
            print("whether a measured folding cost is larger at functionally flagged")
            print("positions at all, across every substitution class in the file.")

    with open(os.path.expanduser(out), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["edge", "direction", "cohort", "window",
                    "n_annotated", "mean_annotated", "cliff_global_annotated",
                    "cliff_matched_annotated", "pathogenic_n_annotated",
                    "n_unannotated", "mean_unannotated", "cliff_global_unannotated",
                    "cliff_matched_unannotated", "pathogenic_n_unannotated"])
        for p in CORE + ARG + EXT:
            r = res.get(p) or {}
            if not r.get('n_ann') and not r.get('n_un'):
                continue
            w.writerow([p[0] + '>' + p[1],
                        'basic' if p in BASICDIR else 'acidic' if (p[1], p[0]) in BASICDIR else '-',
                        a.label, a.window,
                        r.get('n_ann', 0), "%.4f" % r['mean_ann'] if r.get('n_ann') else "",
                        "%.3f" % r['g_ann'] if r.get('n_ann') else "",
                        "%.3f" % r['m_ann'] if r.get('n_ann') else "",
                        r.get('path_ann', ''),
                        r.get('n_un', 0), "%.4f" % r['mean_un'] if r.get('n_un') else "",
                        "%.3f" % r['g_un'] if r.get('n_un') else "",
                        "%.3f" % r['m_un'] if r.get('n_un') else "",
                        r.get('path_un', '')])
    print("\nwrote", os.path.expanduser(out))


if __name__ == "__main__":
    main()
