#!/usr/bin/env python3
"""The substitution set and the single-nucleotide reachability test.

Self-contained. The analysis scripts in this repository import three names from
here that they previously imported from the design method's own modules, which
are not distributed:

    CORE, ARG, EXT, BASICDIR, CHARGE, SCHEMAS   the substitution set
    routes(aa1, aa2)                            single-nucleotide reachability

The twelve core substitutions, the arginine substitutions and the proposed
extensions are all listed below with their codon changes. Nothing here is
taken from the design method's implementation.

routes() is computed from the standard genetic code at import time instead of
being tabulated, so it can be checked by reading twenty lines rather than
trusted. It returns, for an ordered pair of amino acids, every way the first
can become the second by changing one base of one codon, as
(codon1, position, codon2, "X>Y", "ts" or "tv"). The analysis scripts use only
whether that list is empty.

Written 10 October 2026.
"""

BASES = "UCAG"
_AA = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")

CODE = {}
_i = 0
for _b1 in BASES:
    for _b2 in BASES:
        for _b3 in BASES:
            CODE[_b1 + _b2 + _b3] = _AA[_i]
            _i += 1

CODONS = {}
for _c, _a in CODE.items():
    if _a != "*":
        CODONS.setdefault(_a, []).append(_c)

_PURINE = set("AG")


def _kind(a, b):
    """A transition keeps the base within purines or within pyrimidines."""
    return "ts" if (a in _PURINE) == (b in _PURINE) else "tv"


def routes(aa1, aa2):
    """Every single-base route from aa1 to aa2 in the standard genetic code."""
    out = []
    for c1 in CODONS.get(aa1, []):
        for c2 in CODONS.get(aa2, []):
            diff = [i for i in range(3) if c1[i] != c2[i]]
            if len(diff) != 1:
                continue
            i = diff[0]
            out.append((c1, i + 1, c2, "%s>%s" % (c1[i], c2[i]),
                        _kind(c1[i], c2[i])))
    return out


# ------------------------------------------------------- the substitution set
# The twelve core substitutions: the pairs among D, E, H, N, Q and K reachable
# by one change at the first codon position.
CORE = [('D', 'N'), ('N', 'D'), ('E', 'Q'), ('Q', 'E'), ('N', 'H'), ('H', 'N'),
        ('Q', 'K'), ('K', 'Q'), ('D', 'H'), ('H', 'D'), ('E', 'K'), ('K', 'E')]

# The arginine substitutions, reported alongside the twelve.
ARG = [('K', 'R'), ('Q', 'R'), ('R', 'Q'), ('R', 'H'), ('H', 'R')]

# The proposed extensions, reported as outside the first-position rule.
EXT = [('E', 'D'), ('H', 'Y'), ('Q', 'H'), ('D', 'E'), ('H', 'Q'), ('N', 'K'),
       ('K', 'N'), ('Y', 'H')]

# The six core substitutions whose direction raises charge. These are the
# "basic" direction and their reverses the "acidic" one.
BASICDIR = {('E', 'Q'), ('N', 'H'), ('E', 'K'), ('Q', 'K'), ('D', 'N'),
            ('D', 'H')}

# Formal side-chain charge near neutral pH. Histidine is carried at 0.1
# because it is partly protonated in this range.
CHARGE = {'D': -1, 'E': -1, 'K': 1, 'R': 1, 'H': 0.1, 'N': 0, 'Q': 0, 'Y': 0}

# Column names of the two variant cohorts, so a cohort file can be read without
# naming its columns at every call site.
SCHEMAS = [
    {'aa1': 'aa1', 'aa2': 'aa2', 'score': 'score', 'pclass': 'pclass',
     'sector': 'sector', 'snv': 'is_snv'},
    {'aa1': 'a.a.1', 'aa2': 'a.a.2', 'score': 'pathogenicity score',
     'pclass': 'pathogenicity class', 'sector': None, 'snv': 'is_snv'},
]

if __name__ == "__main__":
    print("%d codons, %d amino acids" % (len(CODE), len(CODONS)))
    print("core edges all single-base reachable: %s"
          % all(routes(a, b) for a, b in CORE))
    for a, b in CORE:
        r = routes(a, b)
        print("  %s>%s  %s" % (a, b, ", ".join(
            "%s pos%d %s %s %s" % t for t in r)))
