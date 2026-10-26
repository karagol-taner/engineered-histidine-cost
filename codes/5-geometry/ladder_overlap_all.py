#!/usr/bin/env python3
"""Shape overlap between every pair of residue types.

Applies the measurement in ladder_geometry.py to every type that carries a side
chain and writes the overlap, the change in reach and the change in heavy atoms
for each ordered pair.

    python3 ladder_overlap_all.py <pdbdir> <out.csv>
"""
import csv, glob, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(sys.argv[0])))
import ladder_geometry as LG

THREE = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "ILE", "LEU",
         "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M",
       "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y",
       "VAL": "V"}

if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
    print(__doc__.strip())
    sys.exit(0)
if len(sys.argv) < 3:
    sys.exit("usage: python3 %s <pdbdir> <out.csv>"
             % os.path.basename(sys.argv[0]))
pdbdir, out = sys.argv[1], sys.argv[2]
clouds = {t: [] for t in THREE}
reach = {t: [] for t in THREE}
nfile = 0
for p in sorted(glob.glob(os.path.join(pdbdir, "*.pdb.gz"))):
    atoms = LG.read_pdb(p)
    if not atoms:
        continue
    ch = LG.first_polymer_chain(atoms)
    nfile += 1
    for key, res in LG.residues([a for a in atoms if a[0] == ch]):
        rn = res[0][3]
        if rn not in clouds:
            continue
        loc = LG.sidechain_local(res)
        if loc is None or not len(loc):
            continue
        clouds[rn].append(loc)
        reach[rn].append(float(np.max(np.linalg.norm(loc, axis=1))))
print("structures parsed: %d" % nfile)

occ, n = {}, {}
for t in THREE:
    if not clouds[t]:
        sys.exit("FATAL: no residues of type %s" % t)
    occ[t] = LG.occupancy(np.vstack(clouds[t]))
    n[t] = len(clouds[t])
print("residues measured: %d over %d types" % (sum(n.values()), len(THREE)))
print("  " + "  ".join("%s:%d" % (ONE[t], n[t]) for t in THREE))

rows = [("pair", "overlap", "n_src", "n_dst", "d_reach", "d_atoms")]
med_reach = {t: float(np.median(reach[t])) for t in THREE}
med_atoms = {t: float(np.median([len(c) for c in clouds[t]])) for t in THREE}
for a in THREE:
    for b in THREE:
        if a == b:
            continue
        rows.append((ONE[a] + ">" + ONE[b], "%.4f" % LG.overlap(occ[a], occ[b]),
                     n[a], n[b], "%+.3f" % (med_reach[b] - med_reach[a]),
                     "%+d" % int(med_atoms[b] - med_atoms[a])))
with open(out, "w", newline="") as fh:
    csv.writer(fh).writerows(rows)
print("wrote %s: %d ordered pairs" % (out, len(rows) - 1))
