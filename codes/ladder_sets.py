#!/usr/bin/env python3

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
    return "ts" if (a in _PURINE) == (b in _PURINE) else "tv"


def routes(aa1, aa2):
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


CORE = [('D', 'N'), ('N', 'D'), ('E', 'Q'), ('Q', 'E'), ('N', 'H'), ('H', 'N'),
        ('Q', 'K'), ('K', 'Q'), ('D', 'H'), ('H', 'D'), ('E', 'K'), ('K', 'E')]

ARG = [('K', 'R'), ('Q', 'R'), ('R', 'Q'), ('R', 'H'), ('H', 'R')]

EXT = [('E', 'D'), ('H', 'Y'), ('Q', 'H'), ('D', 'E'), ('H', 'Q'), ('N', 'K'),
       ('K', 'N'), ('Y', 'H')]

BASICDIR = {('E', 'Q'), ('N', 'H'), ('E', 'K'), ('Q', 'K'), ('D', 'N'),
            ('D', 'H')}

CHARGE = {'D': -1, 'E': -1, 'K': 1, 'R': 1, 'H': 0.1, 'N': 0, 'Q': 0, 'Y': 0}

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
