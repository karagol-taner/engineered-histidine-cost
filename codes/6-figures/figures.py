#!/usr/bin/env python3
"""Figures for the engineered-histidine-cost analysis.

Every value is read from the recorded analysis output under data/. Categorical
colour is Okabe-Ito; each series also carries its own marker or line style.

    python3 figures.py
"""
import csv, os, statistics as st, sys

if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
    print(__doc__.strip())
    sys.exit(0)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from matplotlib.colors import to_rgba, LinearSegmentedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.environ.get("LADDER_ROOT") or os.path.normpath(os.path.join(HERE, "..", ".."))
FIG = os.path.join(D, "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 7.5, "axes.labelsize": 7.5, "axes.titlesize": 8.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.edgecolor": "0.35", "grid.color": "0.88", "grid.linewidth": 0.5,
    "grid.linestyle": "-", "savefig.dpi": 600, "figure.dpi": 150,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "axes.unicode_minus": False,
})
INK, MID, LITE = "0.10", "0.45", "0.72"
# Okabe-Ito categorical order.
C1, C2, C3, C4, C5, C6 = "#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9"
SEQ = ("#CCE1F0", "#6BAED6", "#08519C")   # ordinal ramp, light to dark, one hue


def tint(c, a=0.32):
    """A light wash of a categorical hue, for fills larger than a marker."""
    r, g, b, _ = to_rgba(c)
    return (1 - a + a * r, 1 - a + a * g, 1 - a + a * b)


def die(m):
    sys.exit("FATAL: " + m)


def read_core():
    p = os.path.join(D, "data", "ladder_ddg_core.csv")
    if not os.path.exists(p):
        die("missing %s; run ladder_ddg.py first" % p)
    out = {}
    with open(p, newline="") as fh:
        for r in csv.DictReader(fh):
            out[r["edge"]] = r
    if len(out) != 12:
        die("ladder_ddg_core.csv holds %d edges, expected 12" % len(out))
    return out


def parse_block(path, val_i, n_expected=12):
    """Parse the '12 core edges' block of a ladder_edge_cost / ladder_humsavar run."""
    if not os.path.exists(path):
        die("missing %s" % path)
    vals, inblock = {}, False
    for line in open(path):
        if line.startswith("--- 12 core edges"):
            inblock = True
            continue
        if inblock and line.startswith("---"):
            break
        f = line.split()
        if inblock and len(f) > val_i and ">" in f[0]:
            try:
                vals[f[0]] = [float(x) for x in (f[val_i] if isinstance(val_i, int)
                                                 else [f[i] for i in val_i])] \
                    if not isinstance(val_i, int) else float(f[val_i])
            except ValueError:
                pass
    if len(vals) != n_expected:
        die("parsed %d of %d core edges from %s" % (len(vals), n_expected, path))
    return vals


def parse_wilson(path):
    vals = {}
    inblock = False
    for line in open(path):
        if line.startswith("--- 12 core edges"):
            inblock = True
            continue
        if inblock and line.startswith("---"):
            break
        f = line.split()
        if inblock and len(f) > 8 and ">" in f[0]:
            try:
                vals[f[0]] = (float(f[5]), float(f[6]), float(f[8]))
            except ValueError:
                pass
    if len(vals) != 12:
        die("parsed %d of 12 Wilson intervals from %s" % (len(vals), path))
    return vals


def read_strat(window):
    p = os.path.join(D, "data", "ladder_functional_transmembrane_w%d.csv" % window)
    if not os.path.exists(p):
        die("missing %s; run ladder_functional.py --label transmembrane --window %d"
            % (p, window))
    out = {}
    with open(p, newline="") as fh:
        for r in csv.DictReader(fh):
            out[r["edge"]] = r
    return out


def read_fireprot():
    """The per-measurement FireProtDB rows, from the derived table."""
    p = os.path.join(D, "data", "ladder_fireprot_positions.csv")
    if not os.path.exists(p):
        die("missing %s" % p)
    rows = []
    with open(p, newline="") as fh:
        for r in csv.DictReader(fh):
            a = r["accessible_surface_area"]
            rows.append((float(r["ddg_kcal_per_mol"]),
                         float(a) if a else float("nan"),
                         r["functionally_flagged"] == "1"))
    if len(rows) != 3420:
        die("expected 3,420 FireProtDB rows, read %d" % len(rows))
    return rows


CORE = read_core()
TM = parse_block(os.path.join(D, "data", "source_runs", "ladder_edge_cost_tm.out"), 6)
WIL = parse_wilson(os.path.join(D, "data", "source_runs", "ladder_humsavar.out"))
S0, S5 = read_strat(0), read_strat(5)
FP = read_fireprot()

MEGA = {e: float(r["cliff_megascale"]) for e, r in CORE.items()}
AMS = {e: float(r["cliff_alphamissense_soluble"]) for e, r in CORE.items()}
HUM = {e: float(r["humsavar_pathogenic_pct"]) for e, r in CORE.items()}
ORDER = sorted(MEGA, key=lambda e: -MEGA[e])          # by measured folding stability
HIGH = ("D>H", "N>H")


def grid(ax, axis="x"):
    ax.grid(True, axis=axis, zorder=0)
    ax.set_axisbelow(True)


# ------------------------------------------------- measured side-chain geometry
def read_geometry():
    """Numbers from ladder_geometry.py."""
    p = os.path.join(D, "data", "ladder_geometry.csv")
    if not os.path.exists(p):
        die("missing %s; run ladder_geometry.py first" % p)
    out = {}
    with open(p, newline="") as fh:
        for r in csv.DictReader(fh):
            out.setdefault(r["analysis"], {})[r["key"]] = r["value"]
    for need in ("overlap", "d_reach", "reach", "matrix", "predicts_cliff"):
        if need not in out:
            die("ladder_geometry.csv carries no %s block" % need)
    return out


def read_clouds():
    """Frame-aligned side-chain atom clouds, one array per residue type."""
    import numpy as np
    p = os.path.join(D, "data", "ladder_geometry_clouds.npz")
    if not os.path.exists(p):
        die("missing %s; run ladder_geometry.py first" % p)
    z = np.load(p)
    return {k: z[k] for k in z.files}


GEO = read_geometry()
CLOUD = read_clouds()
AA3 = {"D": "ASP", "N": "ASN", "H": "HIS", "E": "GLU", "Q": "GLN",
       "K": "LYS", "Y": "TYR", "R": "ARG"}
# Charge change in units of elementary charge.
DQ = {"D>N": 1, "N>D": 1, "E>Q": 1, "Q>E": 1, "N>H": 0, "H>N": 0,
      "Q>K": 1, "K>Q": 1, "D>H": 1, "H>D": 1, "E>K": 2, "K>E": 2}


# ---------------------------------------------------------------- figure 1
def _view_axes(letters, ref):
    """The 2D view that shows these side chains most fully."""
    import numpy as np
    P = np.vstack([CLOUD["tmpl_" + AA3[L]] for L in letters])
    P = P - P.mean(axis=0)
    _, _, Vt = np.linalg.svd(P, full_matrices=False)
    e1 = Vt[0]
    far = max(P, key=lambda q: float(q @ q))
    if float(e1 @ far) < 0:
        e1 = -e1
    Q = CLOUD["tmpl_" + AA3[ref]]
    Q = Q - Q.mean(axis=0)
    Q = Q - np.outer(Q @ e1, e1)
    _, _, Wt = np.linalg.svd(Q, full_matrices=False)
    e2 = Wt[0]
    e2 = e2 - e1 * float(e2 @ e1)
    e2 = e2 / np.linalg.norm(e2)
    return np.vstack([e1, e2])


VIEW = None


def _draw_sidechain(ax, letter, colour, label=None, lw=1.7, ms=4.0, atoms=False):
    """One observed side chain of this residue type, drawn as sticks."""
    import numpy as np
    t = AA3[letter]
    names = list(CLOUD["tmplnames_" + t])
    xyz = CLOUD["tmpl_" + t]
    pos = {n: VIEW @ p for n, p in zip(names, xyz)}
    ca = VIEW @ np.zeros(3)
    ax.plot([ca[0], pos["CB"][0]], [ca[1], pos["CB"][1]],
            color=colour, linewidth=lw, solid_capstyle="round", zorder=3)
    for u, v in BONDS3[t]:
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                color=colour, linewidth=lw, solid_capstyle="round", zorder=3)
    for n, p in pos.items():
        hetero = n[0] in "NO"
        ax.plot([p[0]], [p[1]], marker="o", markersize=ms if hetero else ms * 0.62,
                markerfacecolor="white" if hetero else colour,
                markeredgecolor=colour, markeredgewidth=0.9 if hetero else 0.0,
                zorder=4)
        if atoms and hetero:
            ax.annotate(n, xy=(p[0], p[1]), xytext=(3.2, 3.2),
                        textcoords="offset points", fontsize=5.3, color=colour)
    return Line2D([], [], color=colour, lw=1.6, label=label or letter)


BONDS3 = {
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


def figure1():
    import numpy as np
    global VIEW

    fig = plt.figure(figsize=(7.1, 6.9))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.96, 0.94, 1.02],
                          hspace=0.46, wspace=0.30,
                          left=0.085, right=0.975, top=0.965, bottom=0.065)
    axa = fig.add_subplot(gs[0, :])
    axb = fig.add_subplot(gs[1, 0])
    axc = fig.add_subplot(gs[1, 1])
    axd = fig.add_subplot(gs[2, 0])
    axe = fig.add_subplot(gs[2, 1])

    # (a) the two closed codon triangles
    axa.set_xlim(0.2, 11.1); axa.set_ylim(1.55, 6.45); axa.axis("off")
    axa.set_aspect("equal", adjustable="box")
    tri = [
        (2.6, [("D", 1.2, 5.0), ("N", 4.0, 5.0), ("H", 2.6, 2.6)],
         "NAY codons, third base U or C",
         {("D", "N"): ("G/A  ts", True), ("D", "H"): ("G/C  tv", False),
          ("N", "H"): ("A/C  tv", False)}),
        (8.0, [("E", 6.6, 5.0), ("K", 9.4, 5.0), ("Q", 8.0, 2.6)],
         "NAR codons, third base A or G",
         {("E", "K"): ("G/A  ts", True), ("E", "Q"): ("G/C  tv", False),
          ("Q", "K"): ("C/A  tv", False)}),
    ]
    for cx, nodes, title, lab in tri:
        pos = {n: (x, y) for n, x, y in nodes}
        axa.text(cx, 6.05, title, ha="center", va="center", fontsize=7.0, color=INK)
        for (u, v), (tx, ts) in lab.items():
            x1, y1 = pos[u]; x2, y2 = pos[v]
            axa.add_patch(FancyArrowPatch(
                (x1, y1), (x2, y2), arrowstyle="<->", mutation_scale=7,
                linewidth=1.3 if ts else 0.9, color=C1 if ts else C2,
                shrinkA=13, shrinkB=13, zorder=2))
            mx, my = (x1 + x2) / 2.0, (y1 + y2) / 2.0
            dx = 0.0 if abs(y1 - y2) < 0.1 else (-0.88 if mx < cx else 0.88)
            dy = 0.28 if abs(y1 - y2) < 0.1 else 0.0
            axa.text(mx + dx, my + dy, tx, ha="center", va="center", fontsize=6.3,
                     color=C1 if ts else C2)
        for n, x, y in nodes:
            axa.add_patch(plt.Circle((x, y), 0.40, facecolor="white",
                                     edgecolor=INK, linewidth=0.8, zorder=3))
            axa.text(x, y, n, ha="center", va="center", fontsize=9.5, zorder=4, color=INK)
    axa.plot([5.3, 5.3], [2.25, 5.6], color=LITE, linewidth=0.5)
    axa.text(5.3, 1.98, "no first-position change crosses between them",
             ha="center", va="top", fontsize=6.1, color=MID)
    axa.legend(handles=[Line2D([], [], color=C1, lw=1.3, label="transition"),
                        Line2D([], [], color=C2, lw=0.9, label="transversion")],
               loc="lower center", bbox_to_anchor=(0.5, -0.17), ncol=2,
               handlelength=1.5, columnspacing=1.8)
    axa.set_title("a", loc="left", fontweight="bold", x=0.0, y=0.99)

    # (b, c) the measured average side chain of each triangle, superposed
    for ax, letters, cols, tag, ref in (
            (axb, "DNH", (C1, C2, C3), "b", "H"),
            (axc, "EQK", (C1, C2, C3), "c", "Q")):
        VIEW = _view_axes(letters, ref)
        hs = []
        for L, c in zip(letters, cols):
            hs.append(_draw_sidechain(ax, L, c, atoms=(L == letters[-1]),
                                      label="%s, reach %.2f A"
                                            % (L, float(GEO["reach"]["median_" + L]))))
        ca = VIEW @ np.zeros(3)
        ax.plot([ca[0]], [ca[1]], marker="o", markersize=4.6, color=INK, zorder=5)
        ax.annotate("CA", xy=(ca[0], ca[1]), xytext=(-5, -10),
                    textcoords="offset points", fontsize=6.2, color=INK)
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_xlabel("angstroms")
        ax.set_ylabel("angstroms")
        ax.legend(handles=hs, loc="upper left", bbox_to_anchor=(-0.02, 1.04),
                  handlelength=1.1, labelspacing=0.22, borderpad=0.2)
        ax.set_title(tag, loc="left", fontweight="bold")
        ax.grid(True, zorder=0); ax.set_axisbelow(True)
    axb.text(0.98, 0.03, "one observed side chain, the medoid of the most\n"
             "populated chi1 well; open circles are N or O",
             transform=axb.transAxes, ha="right", va="bottom", fontsize=5.9, color=MID)

    # (d) the overlap matrix
    letters = list("DNHEQKYR")
    M = np.array([[float(GEO["matrix"][u + v]) for v in letters] for u in letters])
    cmap = LinearSegmentedColormap.from_list("seq", ["#FFFFFF"] + list(SEQ))
    im = axd.imshow(M, cmap=cmap, vmin=0.25, vmax=1.0, zorder=2)
    axd.set_xticks(range(8)); axd.set_xticklabels(letters)
    axd.set_yticks(range(8)); axd.set_yticklabels(letters)
    axd.tick_params(length=0)
    for i in range(8):
        for j in range(8):
            axd.text(j, i, ("%.2f" % M[i, j]).lstrip("0") if M[i, j] < 1 else "1",
                     ha="center", va="center", fontsize=5.6,
                     color="white" if M[i, j] > 0.72 else INK, zorder=3)
    for s in axd.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=axd, orientation="horizontal", fraction=0.050,
                      pad=0.085)
    cb.set_label("shape overlap", fontsize=6.6, labelpad=2)
    cb.ax.tick_params(labelsize=6, width=0.4, length=2)
    cb.outline.set_linewidth(0.4)
    axd.set_title("d", loc="left", fontweight="bold")

    # (e) measured geometry against the measured ordering
    ramp = {0: SEQ[0], 1: SEQ[1], 2: SEQ[2]}
    for e in MEGA:
        axe.scatter([abs(float(GEO["d_reach"][e]))], [MEGA[e]], s=26, marker="o",
                    facecolor=ramp[DQ[e]], edgecolor=INK, linewidth=0.5, zorder=3)
    for e, dx, dy in (("D>H", 9, -2), ("N>H", 9, 2), ("H>D", 9, 1),
                      ("Q>E", 9, -2), ("E>K", -9, 4)):
        axe.annotate(e, xy=(abs(float(GEO["d_reach"][e])), MEGA[e]),
                     xytext=(dx, dy), textcoords="offset points", fontsize=6.4,
                     color=INK, ha="left" if dx > 0 else "right")
    axe.axhline(0, color="0.35", linewidth=0.5, zorder=2)
    axe.set_xlim(-0.10, 1.42); axe.set_ylim(-0.22, 0.46)
    axe.set_xlabel("change in side-chain reach, magnitude (A)")
    axe.set_ylabel("Cliff delta, measured\nfolding stability")
    axe.grid(True, zorder=0); axe.set_axisbelow(True)
    axe.text(0.5, 0.97, "Spearman rho %s over the twelve"
             % GEO["predicts_cliff"]["rho_abs_d_reach"], transform=axe.transAxes,
             ha="center", va="top", fontsize=6.6, color=INK)
    axe.text(0.40, 0.44, "N>H and D>H sit at the same\nplace on this axis",
             transform=axe.transAxes, ha="center", va="center", fontsize=6.1,
             color=MID)
    axe.legend(handles=[plt.Rectangle((0, 0), 1, 1, facecolor=ramp[k], edgecolor=INK,
                                      linewidth=0.4, label=lab)
                        for k, lab in ((0, "charge change ~0"), (1, "1 unit"),
                                       (2, "2 units"))],
               loc="lower left", bbox_to_anchor=(-0.02, -0.02), handlelength=1.1,
               handleheight=0.85, labelspacing=0.22, borderpad=0.2)
    axe.set_title("e", loc="left", fontweight="bold")

    out = os.path.join(FIG, "Figure1_substitution_set.png")
    fig.savefig(out); plt.close(fig)
    print("wrote", out)


# ---------------------------------------------------------------- figure 2

def rank_agreement(name):
    """One rank agreement, read from the deposited significance table."""
    p = os.path.join(D, "data", "ladder_significance.csv")
    with open(p, newline="") as fh:
        for r in csv.reader(fh):
            if len(r) >= 3 and r[0] == "rank agreement" and r[1] == name:
                return r[2].strip()
    die("no rank agreement named %r in %s" % (name, p))


def figure2():
    fig, axs = plt.subplots(2, 2, figsize=(7.1, 5.3),
                            gridspec_kw=dict(hspace=0.62, wspace=0.30))
    axa, axb = axs[0]
    axc, axd = axs[1]
    x = list(range(len(ORDER)))
    series = [("Measured folding stability", MEGA, "o", C1, "-"),
              ("AlphaMissense, soluble", AMS, "s", C2, "--"),
              ("AlphaMissense, transmembrane", TM, "^", C3, ":")]
    for lab, dd, mk, col, ls in series:
        axa.plot(x, [dd[e] for e in ORDER], marker=mk, markersize=3.8, linewidth=0.8,
                 linestyle=ls, color=col, markerfacecolor="white", markeredgecolor=col,
                 markeredgewidth=0.9, label=lab, zorder=3)
    axa.axhline(0, color="0.35", linewidth=0.5, zorder=2)
    axa.set_xticks(x); axa.set_xticklabels(ORDER, rotation=90)
    axa.set_ylabel("Cliff delta")
    axa.set_ylim(-0.33, 0.62)
    grid(axa, "y")
    axa.legend(loc="upper center", bbox_to_anchor=(0.5, 1.30), ncol=1,
               handlelength=2.0, labelspacing=0.22)
    axa.text(0.985, 0.955, "D>H is last in both\nvariant readouts and\neighth in folding stability",
             transform=axa.transAxes, ha="right", va="top", fontsize=6.3, color=INK)
    axa.set_title("a", loc="left", fontweight="bold", y=1.33)

    for i, e in enumerate(ORDER):
        pct, lo, hi = WIL[e]
        axb.plot([i, i], [lo, hi], color=C1, linewidth=0.9, zorder=3)
        axb.plot([i], [pct], marker="o", markersize=3.8, color=C1,
                 markerfacecolor="white", markeredgewidth=1.0, zorder=4)
    axb.axhline(45.5, color="0.35", linewidth=0.6, zorder=2)
    axb.text(0.1, 46.4, "file rate 45.5%", fontsize=6.3, color=MID)
    axb.set_xticks(x); axb.set_xticklabels(ORDER, rotation=90)
    axb.set_ylabel("pathogenic in curated\nhuman variation (%)")
    grid(axb, "y")
    axb.set_title("b", loc="left", fontweight="bold", y=1.33)

    def rk(dd):
        o = sorted(dd, key=lambda e: -dd[e])
        return {e: o.index(e) + 1 for e in dd}
    r_ams, r_tm, r_meg = rk(AMS), rk(TM), rk(MEGA)
    axc.plot([0.5, 12.5], [0.5, 12.5], color=LITE, linewidth=0.5, zorder=1)
    axc.scatter([r_ams[e] for e in ORDER], [r_tm[e] for e in ORDER], s=22, marker="^",
                facecolor="white", edgecolor=C3, linewidth=1.1, zorder=3,
                label="AlphaMissense soluble against\ntransmembrane, rho %s"
                      % rank_agreement("AM soluble vs AM transmembrane"))
    axc.scatter([r_meg[e] for e in ORDER], [r_ams[e] for e in ORDER], s=18, marker="o",
                facecolor=C2, edgecolor="white", linewidth=0.6, zorder=2,
                label="folding stability against\nAlphaMissense soluble, rho %s"
                      % rank_agreement("megascale vs AM soluble"))
    axc.set_xlim(0.5, 12.5); axc.set_ylim(0.5, 12.5)
    axc.set_aspect("equal", adjustable="box")
    axc.set_xlabel("rank in the first readout")
    axc.set_ylabel("rank in the second readout")
    axc.set_xticks([1, 4, 8, 12]); axc.set_yticks([1, 4, 8, 12])
    axc.grid(True, zorder=0); axc.set_axisbelow(True)
    axc.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.22), handlelength=1.0,
               borderpad=0.0, labelspacing=0.5)
    axc.set_title("c", loc="left", fontweight="bold")

    routes = ["N>H", "D>H", "Y>H"]
    yh = {"mega": -0.322, "ams": -0.249, "tm": -0.269}
    vals = {"N>H": (MEGA["N>H"], AMS["N>H"], TM["N>H"]),
            "D>H": (MEGA["D>H"], AMS["D>H"], TM["D>H"]),
            "Y>H": (yh["mega"], yh["ams"], yh["tm"])}
    w = 0.25
    styles = [(C1, "", "Measured folding stability"),
              (C2, "", "AlphaMissense, soluble"),
              (C3, "", "AlphaMissense, transmembrane")]
    for j, r in enumerate(routes):
        for k, (col, hatch, _) in enumerate(styles):
            axd.bar(j + (k - 1) * w, vals[r][k], width=w * 0.9, color=col,
                    edgecolor="white", linewidth=0.6, hatch=hatch, zorder=3)
    axd.axhline(0, color="0.35", linewidth=0.5, zorder=4)
    axd.set_xticks(range(3)); axd.set_xticklabels(routes)
    axd.set_ylabel("Cliff delta")
    axd.set_ylim(-0.40, 0.33)
    grid(axd, "y")
    axd.legend(handles=[plt.Rectangle((0, 0), 1, 1, facecolor=c, hatch=h,
                                      edgecolor=c, linewidth=0.4, label=l)
                        for c, h, l in styles],
               loc="upper left", bbox_to_anchor=(-0.02, -0.14), handlelength=1.2,
               handleheight=0.9, labelspacing=0.3)
    axd.set_title("d", loc="left", fontweight="bold")
    out = os.path.join(FIG, "Figure2_ordering_across_readouts.png")
    fig.savefig(out); plt.close(fig)
    print("wrote", out)


# ---------------------------------------------------------------- figure 3
def _rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def figure3():
    wp_s = _rows(os.path.join(D, "data", "ladder_within_protein_soluble.csv"))
    wp_t = _rows(os.path.join(D, "data", "ladder_within_protein_transmembrane.csv"))
    rb = {r["edge"]: r for r in _rows(os.path.join(D, "data",
                                                   "ladder_robustness_transmembrane.csv"))}

    fig, axs = plt.subplots(2, 2, figsize=(7.1, 5.6),
                            gridspec_kw=dict(hspace=0.56, wspace=0.34))
    axa, axb = axs[0]
    axc, axd = axs[1]

    # (a) every protein's own ordering against the pooled one
    rng = __import__("random").Random(7)
    for rows, lab, face, edge, y0 in (
            (wp_t, "transmembrane, 376 proteins", to_rgba(C1, 0.45),
             to_rgba(C1, 0.85), 1.0),
            (wp_s, "soluble, 24 proteins", "white", C2, 0.0)):
        v = [float(r["spearman_rho_vs_pooled"]) for r in rows]
        ys = [y0 + (rng.random() - 0.5) * 0.42 for _ in v]
        axa.scatter(v, ys, s=9, facecolor=face, edgecolor=edge, linewidth=0.5,
                    zorder=3, label=lab)
        m = st.median(v)
        axa.plot([m, m], [y0 - 0.30, y0 + 0.30], color=INK, linewidth=1.6, zorder=4)
        axa.text(m, y0 + 0.36, "median %+.3f" % m, ha="center", fontsize=6.3, color=INK)
    axa.axvline(0, color="0.35", linewidth=0.6, zorder=2)
    axa.set_yticks([0.0, 1.0])
    axa.set_yticklabels(["soluble", "transmembrane"])
    axa.set_ylim(-0.55, 1.65)
    axa.set_xlim(-0.08, 1.02)
    axa.set_xlabel("Spearman rho against the pooled ordering")
    grid(axa)
    axa.text(0.02, 0.03, "every protein positive", transform=axa.transAxes,
             fontsize=6.3, color=MID)
    axa.set_title("a", loc="left", fontweight="bold")

    # (b) cluster bootstrap over proteins
    es = sorted(rb, key=lambda e: -float(rb[e]["pooled_cliff"]))
    for i, e in enumerate(es):
        lo, hi = float(rb[e]["boot_lo"]), float(rb[e]["boot_hi"])
        c = float(rb[e]["pooled_cliff"])
        axb.plot([lo, hi], [i, i], color=tint(C1, 0.72), linewidth=1.2, zorder=3)
        axb.plot([c], [i], marker="o", markersize=3.8, color=C1,
                 markerfacecolor="white", markeredgewidth=0.9, zorder=4)
    axb.set_yticks(range(len(es))); axb.set_yticklabels(es)
    axb.invert_yaxis()
    axb.axvline(0, color="0.35", linewidth=0.5, zorder=2)
    axb.set_xlabel("Cliff delta, 95% interval over proteins")
    grid(axb)
    axb.text(0.98, 0.04, "D>H last in 100%\nof 1,000 resamples", transform=axb.transAxes,
             ha="right", va="bottom", fontsize=6.3, color=INK)
    axb.set_title("b", loc="left", fontweight="bold")

    # (c) constraint deciles
    show = [("Q>E", "o", C1, "-"), ("N>H", "s", C2, "--"),
            ("H>D", "^", C3, ":"), ("D>H", "D", C5, "-.")]
    xs = list(range(1, 11))
    for e, mk, col, ls in show:
        ys = [float(rb[e]["constraint_decile_%d" % k]) for k in xs]
        axc.plot(xs, ys, marker=mk, markersize=3.6, linewidth=0.9, linestyle=ls,
                 color=col, markerfacecolor="white", markeredgecolor=col,
                 markeredgewidth=0.8, label=e, zorder=3)
    axc.axhline(0, color="0.35", linewidth=0.5, zorder=2)
    axc.set_xticks(xs)
    axc.set_xlabel("decile of positional constraint, 1 = least constrained")
    axc.set_ylabel("matched Cliff delta")
    grid(axc, "y")
    axc.legend(loc="upper center", bbox_to_anchor=(0.5, 1.26), ncol=4, handlelength=1.8,
               columnspacing=1.0)
    axc.set_title("c", loc="left", fontweight="bold")

    # (d) the functional stratification, at an annotation window of zero
    es0 = [e for e in sorted(S0, key=lambda q: -float(S0[q]["cliff_matched_annotated"]))
           if e in MEGA]
    for i, e in enumerate(es0):
        a = float(S0[e]["cliff_matched_annotated"])
        u = float(S0[e]["cliff_matched_unannotated"])
        axd.plot([a, u], [i, i], color=LITE, linewidth=0.8, zorder=2)
        axd.plot([a], [i], marker="o", markersize=3.8, color=C1,
                 markerfacecolor=C1, zorder=3)
        axd.plot([u], [i], marker="o", markersize=3.8, color=C2,
                 markerfacecolor="white", markeredgewidth=0.9, zorder=3)
    axd.set_yticks(range(len(es0))); axd.set_yticklabels(es0)
    axd.invert_yaxis()
    axd.set_xlim(-0.44, 0.70)
    axd.axvline(0, color="0.35", linewidth=0.5, zorder=4)
    axd.set_xlabel("matched Cliff delta")
    grid(axd)
    axd.legend(handles=[Line2D([], [], marker="o", color=C1, markerfacecolor=C1,
                               lw=0, markersize=3.8, label="annotated"),
                        Line2D([], [], marker="o", color=C2, markerfacecolor="white",
                               lw=0, markersize=3.8, label="unannotated")],
               loc="upper center", bbox_to_anchor=(0.5, 1.26), ncol=2, handlelength=0.8)
    axd.text(0.98, 0.04, "D>H last in both strata", transform=axd.transAxes,
             ha="right", fontsize=6.3, color=MID)
    axd.set_title("d", loc="left", fontweight="bold")

    out = os.path.join(FIG, "Figure3_robustness.png")
    fig.savefig(out); plt.close(fig)
    print("wrote", out)


# ---------------------------------------------------------------- figure 4
def figure4():
    v = {}
    with open(os.path.join(D, "data", "ladder_variance.csv"), newline="") as fh:
        for r in csv.DictReader(fh):
            v.setdefault(r["analysis"], {})[r["key"]] = r["value"]
    dz = _rows(os.path.join(D, "data", "ladder_designed.csv"))
    nat = {r["edge"]: r for r in dz if r["class"] == "natural"}
    des = {r["edge"]: r for r in dz if r["class"] == "designed"}

    fig, axs = plt.subplots(2, 2, figsize=(7.1, 5.6),
                            gridspec_kw=dict(hspace=0.62, wspace=0.40))
    axa, axb = axs[0]
    axc, axd = axs[1]

    # (a) where the per-protein effect of a substitution lives
    sets = [("megascale_all", "all 361 classes\nmeasured stability"),
            ("cohort_core", "the twelve\ntransmembrane cohort"),
            ("megascale_core", "the twelve\nmeasured stability")]
    for i, (k, lab) in enumerate(sets):
        sh = float(v[k]["share"])
        axa.barh(i, sh, height=0.52, color=C1, edgecolor="white", linewidth=1.0, zorder=3)
        axa.barh(i, 1 - sh, left=sh, height=0.52, color=tint(C2, 0.45),
                 edgecolor="white", linewidth=1.0, zorder=3)
        axa.text(sh / 2, i, "%.0f%%" % (100 * sh), ha="center", va="center",
                 fontsize=7.2, color="white", zorder=4)
    axa.set_yticks(range(len(sets)))
    axa.set_yticklabels(["%s\n(transfers at R$^2$ %.2f)"
                         % (l, float(v[k + "_transfer"]["r2"]))
                         for k, l in sets])
    axa.invert_yaxis()
    axa.set_xlim(0, 1.0)
    axa.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    axa.set_xticklabels(["0", "25", "50", "75", "100"])
    axa.set_xlabel("share of the per-protein effect variance (%)")
    grid(axa)
    axa.legend(handles=[plt.Rectangle((0, 0), 1, 1, facecolor=C1, edgecolor=C1,
                                      linewidth=0.5, label="the substitution"),
                        plt.Rectangle((0, 0), 1, 1, facecolor=tint(C2, 0.45),
                                      edgecolor=C2, linewidth=0.5,
                                      label="protein by substitution")],
               loc="upper center", bbox_to_anchor=(0.5, 1.30), ncol=2, handlelength=1.3,
               handleheight=0.9)
    axa.set_title("a", loc="left", fontweight="bold", y=1.32)

    # (b) predicting held-out domains from the substitution, the position, or both
    keys = [("headtohead:all substitution classes", "all classes"),
            ("headtohead:the twelve core substitutions", "the twelve")]
    preds = [("r2_sub", "substitution alone", C1), ("r2_pos", "position alone", C2),
             ("r2_both", "both", C3)]
    w = 0.26
    for j, (k, lab) in enumerate(keys):
        for m, (fld, _, col) in enumerate(preds):
            axb.bar(j + (m - 1) * w, float(v[k][fld]), width=w * 0.9, color=col,
                    edgecolor="white", linewidth=0.6, zorder=3)
    axb.set_xticks(range(len(keys)))
    axb.set_xticklabels([l for _, l in keys])
    axb.set_ylabel("R$^2$, held-out domains")
    axb.set_ylim(0, 0.80)
    grid(axb, "y")
    axb.annotate("+0.043", xy=(1 + w, 0.555), xytext=(0, 6), textcoords="offset points",
                 ha="center", fontsize=6.4, color=INK)
    axb.legend(handles=[plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor=c,
                                      linewidth=0.45, label=l)
                        for _, l, c in preds],
               loc="upper center", bbox_to_anchor=(0.5, 1.30), ncol=3, handlelength=1.1,
               handleheight=0.9, columnspacing=0.9)
    axb.set_title("b", loc="left", fontweight="bold", y=1.32)

    # (c) natural against designed domains
    common = [e for e in nat if e in des]
    lab = {"N>H": (9, 4), "D>H": (8, -11), "Y>H": (11, 5)}
    for e in common:
        x, y = float(nat[e]["cliff"]), float(des[e]["cliff"])
        big = e in lab
        axc.scatter([x], [y], s=26 if big else 14, marker="o",
                    facecolor=C1 if big else "white",
                    edgecolor="white" if big else C1,
                    linewidth=0.7 if big else 0.9, zorder=4 if big else 3)
        if big:
            axc.annotate(e, xy=(x, y), xytext=lab[e], textcoords="offset points",
                         fontsize=6.8, color=INK)
    axc.plot([-0.45, 0.60], [-0.45, 0.60], color=LITE, linewidth=0.5, zorder=1)
    axc.axhline(0, color="0.35", linewidth=0.5, zorder=2)
    axc.axvline(0, color="0.35", linewidth=0.5, zorder=2)
    axc.set_xlim(-0.45, 0.45); axc.set_ylim(-0.48, 0.42)
    axc.set_aspect("equal", adjustable="box")
    axc.set_xlabel("Cliff delta, natural domains")
    axc.set_ylabel("Cliff delta, designed domains")
    axc.grid(True, zorder=0); axc.set_axisbelow(True)
    axc.text(0.03, 0.96, "below the diagonal is\nmore costly when designed",
             transform=axc.transAxes, va="top", fontsize=6.2, color=MID)
    axc.set_title("c", loc="left", fontweight="bold")

    # (d) the size of the conservation confound the matched baseline removes
    es0 = [e for e in S0 if e in MEGA]
    for e in es0:
        axd.scatter([float(S0[e]["cliff_global_annotated"])],
                    [float(S0[e]["cliff_matched_annotated"])], s=18, marker="o",
                    facecolor="white", edgecolor=C3, linewidth=1.1, zorder=3)
    axd.annotate("D>H", xy=(float(S0["D>H"]["cliff_global_annotated"]),
                            float(S0["D>H"]["cliff_matched_annotated"])),
                 xytext=(14, -4), textcoords="offset points", fontsize=6.8, color=INK,
                 arrowprops=dict(arrowstyle="-", lw=0.5, color=MID))
    axd.plot([-0.85, 0.70], [-0.85, 0.70], color=LITE, linewidth=0.5, zorder=1)
    axd.axhline(0, color="0.35", linewidth=0.5, zorder=2)
    axd.set_xlim(-0.85, 0.22); axd.set_ylim(-0.45, 0.70)
    axd.set_xlabel("global Cliff delta at annotated positions")
    axd.set_ylabel("matched Cliff delta\nat annotated positions")
    axd.grid(True, zorder=0); axd.set_axisbelow(True)
    axd.text(0.03, 0.95, "distance from the diagonal is the\nconservation confound the "
             "matched\nbaseline removes", transform=axd.transAxes, va="top",
             fontsize=6.2, color=MID)
    axd.set_title("d", loc="left", fontweight="bold")

    out = os.path.join(FIG, "Figure4_scope.png")
    fig.savefig(out); plt.close(fig)
    print("wrote", out)


# ------------------------------------------------------------ figure S1
def figureS1():
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.1, 3.2),
                                   gridspec_kw=dict(wspace=0.36))

    es = [e for e in sorted(S5, key=lambda q: -float(S5[q]["cliff_matched_annotated"]))
          if e in MEGA]
    for i, e in enumerate(es):
        a = float(S5[e]["cliff_matched_annotated"])
        u = float(S5[e]["cliff_matched_unannotated"])
        axa.plot([a, u], [i, i], color=LITE, linewidth=0.8, zorder=2)
        axa.plot([a], [i], marker="o", markersize=3.8, color=C1,
                 markerfacecolor=C1, zorder=3)
        axa.plot([u], [i], marker="o", markersize=3.8, color=C2,
                 markerfacecolor="white", markeredgewidth=0.9, zorder=3)
    axa.set_yticks(range(len(es))); axa.set_yticklabels(es)
    axa.invert_yaxis()
    axa.set_xlim(-0.44, 0.70)
    axa.axvline(0, color="0.35", linewidth=0.5, zorder=4)
    axa.set_xlabel("matched Cliff delta")
    grid(axa)
    axa.legend(handles=[Line2D([], [], marker="o", color=C1, markerfacecolor=C1,
                               lw=0, markersize=3.8, label="annotated"),
                        Line2D([], [], marker="o", color=C2, markerfacecolor="white",
                               lw=0, markersize=3.8, label="unannotated")],
               loc="upper center", bbox_to_anchor=(0.5, 1.22), ncol=2, handlelength=0.8)
    axa.text(0.98, 0.04, "annotation window 5", transform=axa.transAxes,
             ha="right", fontsize=6.3, color=MID)
    axa.set_title("a", loc="left", fontweight="bold")

    fl = [v for v, _, f in FP if f]
    nf = [v for v, _, f in FP if not f]
    bp = axb.boxplot([fl, nf], widths=0.42, showfliers=False, patch_artist=True,
                     medianprops=dict(color=INK, linewidth=1.1),
                     boxprops=dict(linewidth=0.5, edgecolor=INK),
                     whiskerprops=dict(linewidth=0.5, color=MID),
                     capprops=dict(linewidth=0.5, color=MID))
    bp["boxes"][0].set_facecolor(tint(C1, 0.34)); bp["boxes"][0].set_edgecolor(C1)
    bp["boxes"][1].set_facecolor(tint(C2, 0.34)); bp["boxes"][1].set_edgecolor(C2)
    axb.set_xticks([1, 2])
    axb.set_xticklabels(["flagged functional\n(n=%d)" % len(fl),
                         "not flagged\n(n=%d)" % len(nf)])
    axb.set_ylabel("measured folding cost (kcal/mol)")
    axb.set_ylim(-3.8, 11.8)
    axb.axhline(0, color="0.35", linewidth=0.5, zorder=2)
    grid(axb, "y")
    axb.text(0.5, 0.44, "median %+.3f against %+.3f" % (st.median(fl), st.median(nf)),
             transform=axb.transAxes, ha="center", fontsize=6.4, color=MID)
    axb.set_title("b", loc="left", fontweight="bold")
    ins = axb.inset_axes([0.10, 0.70, 0.86, 0.24])
    af = [a for _, a, f in FP if f and a == a]
    an = [a for _, a, f in FP if not f and a == a]
    ins.hist(af, bins=14, density=True, histtype="step", linewidth=0.9, color=C1)
    ins.hist(an, bins=14, density=True, histtype="step", linewidth=0.9,
             color=C2, linestyle="--")
    ins.set_xlabel("accessible surface area", fontsize=5.8, labelpad=1)
    ins.set_yticks([])
    ins.tick_params(labelsize=5.2, width=0.4, length=1.6, pad=1)
    for sp in ("top", "right"):
        ins.spines[sp].set_visible(False)
    ins.text(0.98, 0.95, "flagged (solid), not flagged (dashed)", transform=ins.transAxes,
             ha="right", va="top", fontsize=5.6, color=MID)

    out = os.path.join(FIG, "FigureS1_functional_stratification.png")
    fig.savefig(out); plt.close(fig)
    print("wrote", out)


figure1(); figure2(); figure3(); figure4(); figureS1()
