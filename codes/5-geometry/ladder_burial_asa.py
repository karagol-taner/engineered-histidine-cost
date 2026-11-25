#!/usr/bin/env python3
import argparse, csv, glob, importlib.util, os, sys

import numpy as np
from Bio.PDB.Polypeptide import protein_letters_3to1 as THREE


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


LG = _load(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "ladder_geometry.py"), "ladder_geometry")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdbdir", default="../../data/pdb")
    ap.add_argument("--out", default="../../data/ladder_burial_asa.csv")
    a = ap.parse_args()

    done = set()
    if os.path.exists(a.out):
        with open(a.out, newline="", encoding="utf-8") as fh:
            done = {r["pdb"] for r in csv.DictReader(fh)}
    files = sorted(glob.glob(os.path.join(a.pdbdir, "*.pdb.gz")))
    if not files:
        sys.exit("no structures in %s" % a.pdbdir)
    new = [f for f in files if os.path.basename(f).split(".")[0].upper() not in done]

    first = not os.path.exists(a.out)
    skipped = []
    with open(a.out, "a", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        if first:
            w.writerow(["pdb", "chain", "order", "resname1", "resseq", "asa"])
        for f in new:
            pid = os.path.basename(f).split(".")[0].upper()
            atoms = LG.read_pdb(f)
            if not atoms:
                skipped.append((pid, "no heavy atoms parsed"))
                continue
            ch = LG.first_polymer_chain(atoms)
            atoms = [x for x in atoms if x[0] == ch]
            rs = list(LG.residues(atoms))
            if len(rs) < 20:
                skipped.append((pid, "largest chain under 20 residues"))
                continue
            xyz = np.array([x[6] for x in atoms])
            rad = np.array([LG.VDW[x[5]] for x in atoms])
            per_atom = LG.sasa_per_atom(xyz, rad)
            idx, pos = 0, {}
            for k, ra in rs:
                pos[k] = (idx, idx + len(ra))
                idx += len(ra)
            order = 0
            for k, ra in rs:
                rn = "MET" if ra[0][3] == "MSE" else ra[0][3]
                aa = THREE.get(rn, "")
                if not aa:
                    continue
                i0, i1 = pos[k]
                sc = [j for j in range(i0, i1) if atoms[j][4] not in LG.BACK]
                area = float(per_atom[sc].sum()) if sc else 0.0
                w.writerow([pid, ch, order, aa, k[1], round(area, 2)])
                order += 1
    with open(a.out, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    print("structures %d, cached %d, residues %d, skipped %d"
          % (len(files), len({r["pdb"] for r in rows}), len(rows), len(skipped)))
    for pid, why in skipped:
        print("   skipped %s: %s" % (pid, why))


if __name__ == "__main__":
    main()
