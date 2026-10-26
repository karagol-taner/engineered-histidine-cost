#!/usr/bin/env python3
"""Gene-clustered intervals on the curated-variation pathogenic fractions.

Resamples genes with replacement and takes the percentile interval of the
resampled fractions, for each substitution and for the difference between two.

    python3 ladder_humsavar_clusters.py <humsavar.txt> [out.csv]
"""
import csv, math, os, re, sys
import numpy as np

if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
    print(__doc__.strip())
    sys.exit(0)
if len(sys.argv) < 2:
    sys.exit("usage: python3 %s <humsavar.txt> [out.csv]"
             % os.path.basename(sys.argv[0]))

T31 = {'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E',
       'Gly':'G','His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F',
       'Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V'}
ROW = re.compile(r'^(\S+)\s+(\S+)\s+(VAR_\d+)\s+p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\s+(LP/P|LB/B|US)\b')
CORE = ["Q>E","K>Q","E>Q","Q>K","N>D","N>H","H>N","K>E","D>N","E>K","H>D","D>H"]
BOOT, SEED = 2000, 20261007

path = sys.argv[1]
per_gene = {}
tot_p = tot_b = 0
for line in open(path, errors="replace"):
    m = ROW.match(line)
    if not m:
        continue
    gene, _, _, w3, _, t3, cat = m.groups()
    if cat == "US":
        continue
    w, t = T31.get(w3), T31.get(t3)
    if not w or not t:
        continue
    e = w + ">" + t
    d = per_gene.setdefault(gene, {})
    a, b = d.get(e, (0, 0))
    if cat == "LP/P":
        d[e] = (a + 1, b); tot_p += 1
    else:
        d[e] = (a, b + 1); tot_b += 1

genes = sorted(per_gene)
print("genes: %d, classified variants: %d (%d LP/P, %d LB/B)"
      % (len(genes), tot_p + tot_b, tot_p, tot_b))
print("file-wide pathogenic fraction: %.1f%%\n" % (100.0 * tot_p / (tot_p + tot_b)))


def frac(sel, e):
    p = b = 0
    for g in sel:
        a, c = per_gene[g].get(e, (0, 0))
        p += a; b += c
    return (100.0 * p / (p + b), p + b) if p + b else (float("nan"), 0)


def wilson(p, n, z=1.96):
    if not n:
        return (float("nan"), float("nan"))
    ph = p / 100.0
    d = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / d
    h = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / d
    return (100 * (c - h), 100 * (c + h))


rng = np.random.default_rng(SEED)
idx = [rng.integers(0, len(genes), len(genes)) for _ in range(BOOT)]
samples = [[genes[i] for i in ix] for ix in idx]

rows = [("edge", "n", "pct", "wilson_lo", "wilson_hi", "gene_lo", "gene_hi",
         "n_genes", "wilson_width", "gene_width")]
print("%-6s %7s %7s %8s %-18s %-18s %6s" %
      ("edge", "n", "genes", "pct", "Wilson 95%", "gene-cluster 95%", "ratio"))
dist = {}
for e in CORE:
    pct, n = frac(genes, e)
    ng = sum(1 for g in genes if e in per_gene[g])
    bs = np.array([frac(s, e)[0] for s in samples])
    bs = bs[~np.isnan(bs)]
    dist[e] = bs
    glo, ghi = np.percentile(bs, [2.5, 97.5])
    wlo, whi = wilson(pct, n)
    print("%-6s %7d %7d %7.1f  %6.1f to %-7.1f %6.1f to %-7.1f %6.2f"
          % (e, n, ng, pct, wlo, whi, glo, ghi, (ghi - glo) / (whi - wlo)))
    rows.append((e, n, "%.1f" % pct, "%.1f" % wlo, "%.1f" % whi,
                 "%.1f" % glo, "%.1f" % ghi, ng,
                 "%.1f" % (whi - wlo), "%.1f" % (ghi - glo)))

d = dist["D>H"] - dist["N>H"]
lo, hi = np.percentile(d, [2.5, 97.5])
obs = frac(genes, "D>H")[0] - frac(genes, "N>H")[0]
print("\nD>H minus N>H in pathogenic fraction")
print("  observed %+.1f points, gene-clustered 95%% interval %+.1f to %+.1f"
      % (obs, lo, hi))
print("  excludes zero: %s" % ("yes" if lo > 0 or hi < 0 else "NO"))
print("  bootstrap p, two sided: %.4f"
      % (2 * min((d <= 0).mean(), (d >= 0).mean())))
rows.append(("D>H-N>H", "", "%+.1f" % obs, "", "", "%+.1f" % lo, "%+.1f" % hi,
             "", "", "%.1f" % (hi - lo)))

out = sys.argv[2] if len(sys.argv) > 2 else "../../data/ladder_humsavar_clusters.csv"
with open(out, "w", newline="") as fh:
    csv.writer(fh).writerows(rows)
print("\nwrote %s" % out)
