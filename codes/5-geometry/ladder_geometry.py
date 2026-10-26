#!/usr/bin/env python3
"""Side-chain geometry of the substitution set.

Measures the eight residue types the set moves between in the deposited
structures of the natural cDNA-proteolysis domains. Each side chain is expressed
in a frame built on its own backbone; the overlap of two types is the sum over
0.5 angstrom voxels of the smaller of their two occupancy probabilities.

    python3 ladder_geometry.py --help
"""
import argparse, csv, glob, gzip, math, os, sys
import numpy as np

# Bondi van der Waals radii, angstroms.
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "SE": 1.90}

# The residues the substitution set moves between.
TYPES = ["ASP", "ASN", "HIS", "GLU", "GLN", "LYS", "TYR", "ARG"]
ONE = {"ASP": "D", "ASN": "N", "HIS": "H", "GLU": "E", "GLN": "Q",
       "LYS": "K", "TYR": "Y", "ARG": "R"}
BACK = {"N", "CA", "C", "O", "OXT"}

# The substitution set, as source and target one-letter codes.
CORE = ["E>Q", "K>Q", "Q>K", "Q>E", "N>H", "D>N", "N>D", "D>H",
        "E>K", "K>E", "H>N", "H>D"]
EXTRA = ["Y>H", "K>R"]

VOXEL = 0.5          # angstroms, the edge of an occupancy voxel
PROBE = 1.4          # angstroms, the solvent probe
SPHERE_N = 194       # points per atom for Shrake and Rupley

# Covalent connectivity of each side chain.
BONDS = {
    "ASP": [("CB", "CG"), ("CG", "OD1"), ("CG", "OD2")],
    "ASN": [("CB", "CG"), ("CG", "OD1"), ("CG", "ND2")],
    "HIS": [("CB", "CG"), ("CG", "ND1"), ("ND1", "CE1"), ("CE1", "NE2"),
            ("NE2", "CD2"), ("CD2", "CG")],
    "GLU": [("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "OE2")],
    "GLN": [("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "NE2")],
    "LYS": [("CB", "CG"), ("CG", "CD"), ("CD", "CE"), ("CE", "NZ")],
    "TYR": [("CB", "CG"), ("CG", "CD1"), ("CD1", "CE1"), ("CE1", "CZ"),
            ("CZ", "CE2"), ("CE2", "CD2"), ("CD2", "CG"), ("CZ", "OH")],
    "ARG": [("CB", "CG"), ("CG", "CD"), ("CD", "NE"), ("NE", "CZ"),
            ("CZ", "NH1"), ("CZ", "NH2")],
}


def die(m):
    sys.exit("FATAL: " + m)


# ----------------------------------------------------------------- parsing
def read_pdb(path):
    """First model, protein heavy atoms, highest-occupancy altloc only."""
    out, seen_model = [], False
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", errors="replace") as fh:
        for line in fh:
            rec = line[:6]
            if rec == "MODEL ":
                if seen_model:
                    break
                seen_model = True
                continue
            if rec == "ENDMDL":
                break
            if rec not in ("ATOM  ", "HETATM"):
                continue
            alt = line[16]
            if alt not in (" ", "A"):
                continue
            name = line[12:16].strip()
            if name.startswith("H") or name[:1].isdigit():
                continue
            resname = line[17:20].strip()
            if rec == "HETATM" and resname != "MSE":
                continue
            el = (line[76:78].strip() or name[0]).upper()
            if el == "H" or el not in VDW:
                continue
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            out.append((line[21], line[22:26].strip(), line[26],
                        resname, name, el, xyz))
    return out


def residues(atoms):
    """Group parsed atoms into residues, keeping deposition order."""
    cur, key = [], None
    for a in atoms:
        k = (a[0], a[1], a[2])
        if k != key:
            if cur:
                yield key, cur
            cur, key = [], k
        cur.append(a)
    if cur:
        yield key, cur


def first_polymer_chain(atoms):
    """The chain carrying the most residues, which is the deposited domain."""
    n = {}
    for k, _ in residues(atoms):
        n[k[0]] = n.get(k[0], 0) + 1
    if not n:
        return None
    return max(n, key=lambda c: (n[c], c))


# -------------------------------------------------------------- the frame
def frame_of(res_atoms):
    """Orthonormal frame on a residue's own backbone, or None if incomplete."""
    p = {a[4]: np.array(a[6]) for a in res_atoms}
    if not {"N", "CA", "CB"} <= set(p):
        return None
    ca = p["CA"]
    e1 = p["CB"] - ca
    n1 = np.linalg.norm(e1)
    if n1 < 1e-6:
        return None
    e1 = e1 / n1
    v = p["N"] - ca
    v = v - e1 * float(v @ e1)
    n2 = np.linalg.norm(v)
    if n2 < 1e-6:
        return None
    e2 = v / n2
    e3 = np.cross(e1, e2)
    return ca, np.vstack([e1, e2, e3])


def sidechain_local(res_atoms):
    """Side-chain heavy atoms from CB outward, in the residue's own frame."""
    fr = frame_of(res_atoms)
    if fr is None:
        return None
    ca, R = fr
    pts = [R @ (np.array(a[6]) - ca) for a in res_atoms if a[4] not in BACK]
    return np.array(pts) if pts else None


def sidechain_named(res_atoms):
    """The same atoms, keyed by name, for the rotamer average."""
    fr = frame_of(res_atoms)
    if fr is None:
        return None
    ca, R = fr
    return {a[4]: R @ (np.array(a[6]) - ca)
            for a in res_atoms if a[4] not in BACK}


def dihedral(p0, p1, p2, p3):
    """Signed torsion in degrees."""
    b0, b1, b2 = p1 - p0, p2 - p1, p3 - p2
    n1, n2 = np.cross(b0, b1), np.cross(b1, b2)
    m = np.cross(n1, b1 / np.linalg.norm(b1))
    return math.degrees(math.atan2(float(m @ n2), float(n1 @ n2)))


def chi1_well(res_atoms):
    """Which of the three staggered chi1 wells this side chain sits in."""
    p = {a[4]: np.array(a[6]) for a in res_atoms}
    if not {"N", "CA", "CB", "CG"} <= set(p):
        return None
    c = dihedral(p["N"], p["CA"], p["CB"], p["CG"])
    if -120 <= c < 0:
        return "g-"
    if 0 <= c < 120:
        return "g+"
    return "t"


# -------------------------------------------- solvent accessible surface area
def sphere_points(n):
    """Near-uniform points on the unit sphere, by the golden spiral."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = math.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)


SPH = sphere_points(SPHERE_N)


def sasa_per_atom(xyz, radii):
    """Shrake and Rupley, vectorised over atoms through a neighbour tree."""
    from scipy.spatial import cKDTree
    r = radii + PROBE
    tree = cKDTree(xyz)
    rmax = float(r.max())
    out = np.zeros(len(xyz))
    neigh = tree.query_ball_point(xyz, rmax + rmax)
    for i in range(len(xyz)):
        js = [j for j in neigh[i] if j != i]
        pts = xyz[i] + r[i] * SPH
        if js:
            d = np.linalg.norm(pts[:, None, :] - xyz[js][None, :, :], axis=2)
            free = np.all(d >= r[js][None, :], axis=1)
        else:
            free = np.ones(len(pts), dtype=bool)
        out[i] = 4 * math.pi * r[i] ** 2 * free.mean()
    return out


# ----------------------------------------------------------------- overlap
def occupancy(cloud, voxel=VOXEL):
    """Voxel probability of a pooled atom cloud."""
    k = np.floor(cloud / voxel).astype(np.int64)
    uniq, cnt = np.unique(k, axis=0, return_counts=True)
    return {tuple(u): c / cnt.sum() for u, c in zip(uniq, cnt)}


def overlap(a, b):
    """Sum over voxels of the smaller of two probabilities."""
    return sum(min(v, b[k]) for k, v in a.items() if k in b)


# -------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdbdir", required=True)
    ap.add_argument("--core", required=True,
                help="per-substitution Cliff delta table; pass "
                     "data/ladder_core_fullprec.csv, the full-precision "
                     "values, not the rounded "
                     "data/ladder_ddg_core.csv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--clouds", required=True)
    ap.add_argument("--min-per-type", type=int, default=200)
    a = ap.parse_args()

    files = sorted(glob.glob(os.path.join(a.pdbdir, "*.pdb.gz"))) + \
            sorted(glob.glob(os.path.join(a.pdbdir, "*.pdb")))
    if not files:
        die("no structures in %s" % a.pdbdir)

    print("ladder_geometry.py\n")
    print("structures found: %d" % len(files))

    clouds = {t: [] for t in TYPES}
    reach = {t: [] for t in TYPES}
    sasa = {t: [] for t in TYPES}
    nres = {t: 0 for t in TYPES}
    named = {t: {"g-": [], "t": [], "g+": []} for t in TYPES}
    used, skipped = [], []

    for f in files:
        pid = os.path.basename(f).split(".")[0].upper()
        atoms = read_pdb(f)
        if not atoms:
            skipped.append((pid, "no heavy atoms parsed"))
            continue
        ch = first_polymer_chain(atoms)
        atoms = [x for x in atoms if x[0] == ch]
        rs = list(residues(atoms))
        if len(rs) < 20:
            skipped.append((pid, "largest chain under 20 residues"))
            continue

        xyz = np.array([x[6] for x in atoms])
        rad = np.array([VDW[x[5]] for x in atoms])
        per_atom = sasa_per_atom(xyz, rad)
        idx, pos = 0, {}
        for k, ra in rs:
            pos[k] = (idx, idx + len(ra))
            idx += len(ra)

        n_here = 0
        for k, ra in rs:
            rn = ra[0][3]
            if rn == "MSE":
                rn = "MET"
            if rn not in TYPES:
                continue
            loc = sidechain_local(ra)
            if loc is None or len(loc) < 2:
                continue
            i0, i1 = pos[k]
            sc = [j for j in range(i0, i1) if atoms[j][4] not in BACK]
            clouds[rn].append(loc)
            reach[rn].append(float(np.linalg.norm(loc, axis=1).max()))
            sasa[rn].append(float(per_atom[sc].sum()))
            w = chi1_well(ra)
            nm = sidechain_named(ra)
            if w and nm and all(x in nm for x, _ in BONDS[rn]) \
                    and all(y in nm for _, y in BONDS[rn]):
                named[rn][w].append(nm)
            nres[rn] += 1
            n_here += 1
        used.append((pid, ch, len(rs), n_here))

    print("structures used: %d, skipped: %d" % (len(used), len(skipped)))
    for pid, why in skipped:
        print("  skipped %s: %s" % (pid, why))

    print("\nresidues measured per type")
    print("%-5s %-6s %8s %10s %10s %10s" %
          ("type", "letter", "n", "atoms", "reach A", "sc SASA A2"))
    for t in TYPES:
        if nres[t] < a.min_per_type:
            die("only %d %s residues, under the floor of %d; the overlap would "
                "be an estimate of too few rotamers to report"
                % (nres[t], t, a.min_per_type))
        na = int(np.median([len(c) for c in clouds[t]]))
        print("%-5s %-6s %8d %10d %10.2f %10.1f"
              % (t, ONE[t], nres[t], na, float(np.median(reach[t])),
                 float(np.median(sasa[t]))))

    pooled = {t: np.vstack(clouds[t]) for t in TYPES}
    occ = {t: occupancy(pooled[t]) for t in TYPES}

    print("\nshape overlap between types, in a common backbone frame")
    rows = [("analysis", "key", "value")]
    for t in TYPES:
        rows.append(("counts", "n_" + ONE[t], nres[t]))
        rows.append(("counts", "atoms_" + ONE[t],
                     int(np.median([len(c) for c in clouds[t]]))))
        rows.append(("reach", "median_" + ONE[t], "%.3f" % float(np.median(reach[t]))))
        rows.append(("sasa", "median_" + ONE[t], "%.2f" % float(np.median(sasa[t]))))

    ov = {}
    print("%-7s %9s %12s %12s" % ("subst", "overlap", "d reach A", "d atoms"))
    for e in CORE + EXTRA:
        s, t = e.split(">")
        rs_, rt_ = [k for k in TYPES if ONE[k] == s][0], [k for k in TYPES if ONE[k] == t][0]
        o = overlap(occ[rs_], occ[rt_])
        dr = float(np.median(reach[rt_])) - float(np.median(reach[rs_]))
        da = int(np.median([len(c) for c in clouds[rt_]])) - \
             int(np.median([len(c) for c in clouds[rs_]]))
        ov[e] = (o, dr, da)
        print("%-7s %9.3f %12.2f %12d" % (e, o, dr, da))
        rows.append(("overlap", e, "%.3f" % o))
        rows.append(("d_reach", e, "%+.3f" % dr))
        rows.append(("d_atoms", e, "%+d" % da))

    # ---- does measured geometry predict the measured ordering?
    cliff = {}
    with open(a.core, newline="") as fh:
        for r in csv.DictReader(fh):
            k = r.get("edge") or r.get("substitution")
            if k:
                for f in ("cliff_megascale", "cliff", "megascale_cliff"):
                    if f in r and r[f] not in ("", None):
                        cliff[k.strip()] = float(r[f])
                        break
    have = [e for e in CORE if e in cliff]
    if len(have) != 12:
        die("expected the Cliff delta of all twelve substitutions in %s, found "
            "%d: %s" % (a.core, len(have), sorted(have)))

    def spearman(x, y):
        def rank(v):
            o = sorted(range(len(v)), key=lambda i: v[i])
            r = [0.0] * len(v)
            i = 0
            while i < len(o):
                j = i
                while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
                    j += 1
                for k in range(i, j + 1):
                    r[o[k]] = (i + j) / 2.0 + 1
                i = j + 1
            return r
        rx, ry = rank(x), rank(y)
        mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
        num = sum((p - mx) * (q - my) for p, q in zip(rx, ry))
        den = math.sqrt(sum((p - mx) ** 2 for p in rx) *
                        sum((q - my) ** 2 for q in ry))
        return num / den if den else float("nan")

    y = [cliff[e] for e in have]
    r_ov = spearman([ov[e][0] for e in have], y)
    r_dr = spearman([abs(ov[e][1]) for e in have], y)
    r_da = spearman([abs(ov[e][2]) for e in have], y)
    print("\nmeasured geometry against the measured ordering, over the twelve")
    print("  shape overlap                Spearman rho %+.3f" % r_ov)
    print("  change in reach, magnitude   Spearman rho %+.3f" % r_dr)
    print("  change in atom count, magn.  Spearman rho %+.3f" % r_da)
    rows += [("predicts_cliff", "rho_overlap", "%.3f" % r_ov),
             ("predicts_cliff", "rho_abs_d_reach", "%.3f" % r_dr),
             ("predicts_cliff", "rho_abs_d_atoms", "%.3f" % r_da),
             ("predicts_cliff", "n", len(have))]

    # ---- is the overlap an artefact of unequal sampling or of the voxel size?
    rng = np.random.default_rng(7)
    nmin = min(nres[t] for t in TYPES)
    reps = {e: [] for e in CORE + EXTRA}
    for _ in range(20):
        sub = {}
        for t in TYPES:
            pick = rng.choice(nres[t], size=nmin, replace=False)
            sub[t] = occupancy(np.vstack([clouds[t][i] for i in pick]))
        for e in CORE + EXTRA:
            s, t = e.split(">")
            rs_ = [k for k in TYPES if ONE[k] == s][0]
            rt_ = [k for k in TYPES if ONE[k] == t][0]
            reps[e].append(overlap(sub[rs_], sub[rt_]))
    print("\nequal-sample check, every type cut to %d residues, 20 repeats" % nmin)
    print("%-7s %9s %9s %9s" % ("subst", "full", "mean", "sd"))
    for e in CORE + EXTRA:
        m, s = float(np.mean(reps[e])), float(np.std(reps[e]))
        print("%-7s %9.3f %9.3f %9.3f" % (e, ov[e][0], m, s))
        rows.append(("overlap_equal_n", e, "%.3f" % m))
        rows.append(("overlap_equal_n_sd", e, "%.3f" % s))
    rows.append(("overlap_equal_n", "n_per_type", nmin))

    print("\nvoxel sensitivity")
    print("%-8s %12s %14s" % ("voxel A", "rho overlap", "rho |d reach|"))
    for vx in (0.4, 0.5, 0.6):
        oc = {t: occupancy(pooled[t], vx) for t in TYPES}
        ovx = {}
        for e in have:
            s, t = e.split(">")
            rs_ = [k for k in TYPES if ONE[k] == s][0]
            rt_ = [k for k in TYPES if ONE[k] == t][0]
            ovx[e] = overlap(oc[rs_], oc[rt_])
        r1 = spearman([ovx[e] for e in have], y)
        print("%-8.1f %12.3f %14.3f" % (vx, r1, r_dr))
        rows.append(("voxel_sensitivity", "rho_overlap_%.1f" % vx, "%.3f" % r1))

    # ---- the full matrix, for the figure
    print("\nfull overlap matrix")
    print("      " + "".join("%7s" % ONE[t] for t in TYPES))
    for u in TYPES:
        line = "%-6s" % ONE[u]
        for v in TYPES:
            o = 1.0 if u == v else overlap(occ[u], occ[v])
            line += "%7.3f" % o
            rows.append(("matrix", ONE[u] + ONE[v], "%.3f" % o))
        print(line)

    # ---- a representative observed side chain, for the figure
    # modal chi1 well, medoid side chain
    print("\nrepresentative side chain, medoid of the most populated chi1 well")
    print("%-5s %-6s %8s %8s %9s   %s"
          % ("type", "well", "n in", "share", "rmsd A", "atoms"))
    tmpl = {}
    for t in TYPES:
        w = max(named[t], key=lambda k: len(named[t][k]))
        pool = [d for d in named[t][w]]
        if len(pool) < 40:
            die("only %d %s side chains in the %s well, too few to represent"
                % (len(pool), t, w))
        names = sorted({k for d in pool for k in d})
        pool = [d for d in pool if all(k in d for k in names)]
        X = np.array([[d[k] for k in names] for d in pool])      # n, atoms, 3
        dif = X[:, None, :, :] - X[None, :, :, :]
        rms = np.sqrt((dif ** 2).sum(axis=3).mean(axis=2))       # n, n
        i = int(rms.sum(axis=1).argmin())
        share = len(named[t][w]) / sum(len(v) for v in named[t].values())
        print("%-5s %-6s %8d %8.2f %9.2f   %s"
              % (t, w, len(pool), share, float(rms[i].mean()), " ".join(names)))
        tmpl[t] = (w, share, names, X[i])
        rows.append(("rotamer", "well_" + ONE[t], w))
        rows.append(("rotamer", "share_" + ONE[t], "%.3f" % share))
        rows.append(("rotamer", "n_" + ONE[t], len(pool)))
        rows.append(("rotamer", "medoid_rmsd_" + ONE[t], "%.3f" % float(rms[i].mean())))

    rows += [("provenance", "structures_used", len(used)),
             ("provenance", "structures_skipped", len(skipped)),
             ("provenance", "voxel_angstrom", "%.2f" % VOXEL),
             ("provenance", "probe_angstrom", "%.2f" % PROBE),
             ("provenance", "sphere_points", SPHERE_N)]

    with open(a.out, "w", newline="") as fh:
        csv.writer(fh).writerows(rows)
    save = {t: pooled[t] for t in TYPES}
    for t in TYPES:
        w, share, names, xyz = tmpl[t]
        save["tmpl_" + t] = xyz
        save["tmplnames_" + t] = np.array(names)
    np.savez_compressed(a.clouds, **save)
    with open(os.path.splitext(a.out)[0] + "_structures.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["pdb", "chain", "residues", "measured_residues"])
        w.writerows(used)
    print("\nwrote %s, %s and %s"
          % (a.out, a.clouds, os.path.splitext(a.out)[0] + "_structures.csv"))


if __name__ == "__main__":
    main()
