#!/usr/bin/env python3
import bisect, collections, csv, glob, importlib.util, os, re, sys

MUT = re.compile(r"^([A-Z])(\d+)([A-Z])$")
AA = set("ACDEFGHIKLMNPQRSTVWY")


def load(path, name):
    if not os.path.exists(path):
        sys.exit("required module not found: %s" % path)
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


def read_megascale(d):
    out = []
    for p in sorted(glob.glob(os.path.join(d, "*.csv"))):
        with open(p, newline="", encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader(fh):
                m = MUT.match((r.get("mut_type") or "").strip())
                if not m:
                    continue
                a1, pos, a2 = m.group(1), int(m.group(2)), m.group(3)
                if a1 not in AA or a2 not in AA or a1 == a2:
                    continue
                try:
                    v = float(r.get("ddG_ML"))
                except (TypeError, ValueError):
                    continue
                name = (r.get("WT_name") or "").strip()
                name = name[:-4] if name.lower().endswith(".pdb") else name
                out.append((name.upper(), pos, a1, a2, v, (r.get("aa_seq") or "").strip()))
    return out


def chain_sequences(path):
    seqs = collections.defaultdict(list)
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            seqs[(r["pdb"], r["chain"])].append((int(r["order"]), r["resname1"], float(r["asa"])))
    out = {}
    for k, v in seqs.items():
        v.sort()
        out[k] = ("".join(a for _, a, _ in v), [s for _, _, s in v])
    return out


def best_offset(dom_seq, chain_seq):
    if not dom_seq or len(dom_seq) > len(chain_seq):
        return None, 0.0
    best, bid = None, 0.0
    for o in range(len(chain_seq) - len(dom_seq) + 1):
        hit = sum(1 for i, c in enumerate(dom_seq) if chain_seq[o + i] == c)
        idn = hit / len(dom_seq)
        if idn > bid:
            best, bid = o, idn
    return best, bid


def rank_edges(rows, snv_ok, CORE):
    base = [v for _, _, a1, a2, v in rows if snv_ok((a1, a2))]
    if len(base) < 200:
        return None, 0
    cliff, nb = cliff_factory(base)
    per = collections.defaultdict(list)
    for _, _, a1, a2, v in rows:
        per[(a1, a2)].append(v)
    out = {}
    for p in CORE:
        if len(per.get(p, [])) >= 20:
            out[p] = (cliff(per[p]), len(per[p]))
    return out, nb
