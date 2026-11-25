# Directional cost of an engineered histidine in near-isosteric charge-altering substitutions

**Taner Karagöl & Alper Karagöl**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Abstract

Proteins working in the bloodstream may have to lose activity in a tumour or an endosome, a unit lower in pH. Histidine, the one amino acid that takes on charge across the physiological range, is introduced to confer that response, and the residue it replaces is chosen by screening. We ranked the twelve near-isosteric charge-altering substitutions reachable by one first-position codon change against four readouts: 271,231 folding-stability values, predicted pathogenicity over cohorts of 306,622 soluble and 5,846,493 transmembrane substitutions, and 73,296 classified human variants. Introducing a histidine is not the reverse of removing one. The two routes into one carry the difference: D>H lies at +0.150 and H>D at -0.117, N>H at +0.266 and H>N at +0.060, gaps of 0.267 and 0.205 whose intervals exclude zero. The other four differ by 0.039, none resolved. BLOSUM62 predicts the ordering at rank correlation +0.669 and the source residue's general tolerance at +0.126. Removing the source effect leaves the asymmetry larger. BLOSUM62, the Grantham distance and shape overlap are symmetric by construction, so each scores a directional gap as zero. The sign of the asymmetry is phenotype-specific. Introducing a histidine on an aspartate is the more damaging direction under a variant readout and removing one the more destabilising under proteolysis, so a design must name the phenotype it selects on. The ordering holds across all four readouts, at Spearman rho +0.622 to +0.783 between the stability data and the variant readouts, and higher among the variant readouts, which share a predictor.

**Keywords:** histidine scanning, pH-responsive protein design, substitution asymmetry, substitution tolerance, variant effect prediction, pI-code

## Repository layout

```
engineered-histidine-cost/
  codes/                             2 files
  codes/1-cost/                      1 file
  codes/2-stratification/            2 files
  codes/3-robustness/                2 files
  codes/4-variance/                  4 files
  codes/5-geometry/                  6 files
  codes/6-figures/                   1 file
  codes/7-significance/              1 file
  data/                              29 files
  data/source_runs/                  16 files
  figures/                           5 files
  LICENSE
  README.md
  Supplementary_Information.pdf
```

* **`codes/1-cost/`**: `ladder_ddg.py` reads the three stability sets, the two AlphaMissense cohorts and the humsavar variant file, maps every stability column to a common sign, and writes the per-substitution Cliff delta for each readout.
* **`codes/2-stratification/`**: `ladder_humsavar_clusters.py` gives the curated-variation fractions a gene-clustered interval. `ladder_functional.py` splits a cohort into positions UniProt annotates as an active site, a binding site or a site and the rest, and scores each substitution against the baseline computed inside its own stratum.
* **`codes/3-robustness/`**: `ladder_robustness.py` ranks the substitutions inside each protein separately, resamples proteins rather than positions for the intervals, and stratifies by a leave-one-out constraint score. `ladder_designed.py` compares designed against natural domains.
* **`codes/4-variance/`**: `ladder_variance.py` separates the per-protein effect into a between-substitution and a protein-by-substitution component and predicts held-out domains. `ladder_extra.py` runs the descriptor regressions, the alternative baselines, the topology comparison and the transition counts. `ladder_source_null.py` puts the source-residue null and BLOSUM62 against the ordering. `ladder_direction.py` takes each reciprocal pair forward against reverse and sets an interval on the gap by a cluster bootstrap over domains.
* **`codes/5-geometry/`**: `ladder_geometry.py` measures the eight residue types in the deposited structures of the natural cDNA-proteolysis domains and writes the shape overlaps, reaches and representative side chains. `ladder_overlap_all.py` widens the same measurement to every side-chain type and `ladder_overlap_stability.py` tests the overlap against measured folding cost over every substitution class. `fetch_structures.py` obtains the depositions by identifier, `ladder_burial_asa.py` writes the side-chain solvent accessibility of every position in them, and `ladder_burial.py` ranks the substitutions separately among buried and among exposed positions over a sweep of accessibility thresholds.
* **`codes/6-figures/`**: `figures.py` rebuilds all five figures from `data/` alone, so they can be checked without repeating the analysis.
* **`codes/7-significance/`**: `ladder_significance.py` tests each substitution against its baseline on the cluster bootstrap and each rank agreement by permutation, with Benjamini-Hochberg applied across the twelve and across the eight comparisons.
* **`data/`**: Every derived table the paper draws a number from, with the printed log of each run in `data/source_runs/`.
* **`figures/`**: The published figures at 600 dpi.

`codes/ladder_sets.py` carries the substitution set and the single-nucleotide reachability test, and `codes/ladder_megascale.py` holds the one reader and scoring rule the cDNA-proteolysis scripts share. Reachability is computed from the standard genetic code at import rather than tabulated, so it can be read rather than trusted. `Supplementary_Information.pdf` holds the supplementary tables and figure.

## Dataset summary

| File | Size | Contents |
| --- | --- | --- |
| `ladder_ddg_core.csv` | 12 rows | Cliff delta of each substitution in every readout, at the precision the tables print. |
| `ladder_core_fullprec.csv` | 12 rows | The same values at full precision, which the reciprocal-pair gaps are taken from. |
| `ladder_source_null.csv` | 12 rows | The source-residue null and BLOSUM62 against the measured ordering. |
| `ladder_direction.csv` | 10 rows | Each reciprocal pair forward against reverse, with a cluster-bootstrap interval on the gap. |
| `ladder_class_costs.csv` | 361 rows | Every substitution class with at least twenty records, which the source-residue null is built from. |
| `ladder_significance.csv` | 39 rows | Each substitution against its baseline and each rank agreement, with Benjamini-Hochberg correction. |
| `ladder_functional.csv` | 25 rows | Annotated against unannotated positions in the soluble cohort. |
| `ladder_functional_transmembrane_w0.csv` | 25 rows | The same stratification in the transmembrane cohort, annotation window zero. |
| `ladder_functional_transmembrane_w5.csv` | 25 rows | The same, annotation window five residues. |
| `ladder_robustness_soluble.csv` | 12 rows | Per-protein ranking and the cluster bootstrap, soluble cohort. |
| `ladder_robustness_transmembrane.csv` | 12 rows | The same for the transmembrane cohort. |
| `ladder_within_protein_soluble.csv` | 24 rows | Each qualifying protein's own ordering against the pooled one, soluble cohort. |
| `ladder_within_protein_transmembrane.csv` | 376 rows | The same for the transmembrane cohort. |
| `ladder_percell_soluble.csv` | 288 rows | The per-protein effect of each substitution, soluble cohort. |
| `ladder_percell_transmembrane.csv` | 4,512 rows | The same for the transmembrane cohort. |
| `ladder_variance.csv` | 51 rows | Variance of the per-protein effect and the held-out prediction. |
| `ladder_designed.csv` | 32 rows | Designed against natural domains in the cDNA-proteolysis set. |
| `ladder_extra.csv` | 10 rows | Descriptor regressions, alternative baselines, topology and transition counts. |
| `ladder_humsavar_clusters.csv` | 13 rows | Curated-variation pathogenic fractions with gene-clustered intervals. |
| `ladder_geometry.csv` | 211 rows | Measured side-chain geometry of the eight residue types. |
| `ladder_geometry_structures.csv` | 188 rows | The deposited identifiers the geometry is measured on. |
| `ladder_geometry_clouds.npz` | 631 KB | Frame-aligned side chains and their templates. |
| `ladder_overlap_all_pairs.csv` | 342 rows | Shape overlap between every pair of residue types. |
| `ladder_overlap_stability.csv` | 10 rows | Shape overlap against measured folding cost. |
| `ladder_burial_asa.csv` | 11,694 rows | Side-chain solvent accessibility at every position of the natural domains. |
| `ladder_burial.csv` | 20 rows | The ordering inside buried and exposed strata over a sweep of thresholds. |
| `ladder_fireprot_positions.csv` | 3,420 rows | Per-measurement folding cost, accessible surface area and functional flag. |
| `cohort_soluble_66.csv` | 66 rows | UniProt accessions of the soluble AlphaMissense cohort, with the protein index. |
| `cohort_transmembrane_504.csv` | 504 rows | UniProt accessions of the transmembrane AlphaMissense cohort. |

The derived tables are the output of the runs in `data/source_runs/`, so every reported number and every figure can be checked without repeating them. Figures 1 to 4 and Figure S1 are rebuilt from these tables alone.

## Reproducing the analysis

Each script runs from its own folder and finds the repository from its path, or from `LADDER_ROOT` when that is set. The scripts that take options show their inputs under `--help`, and their defaults expect the primary files beside the derived tables, as `data/megascale/`, `data/fireprot/`, `data/s669/`, `data/humsavar.txt`, `data/alphamissense_soluble.csv.gz` and `data/alphamissense_transmembrane/`.

```
cd codes/6-figures
python3 figures.py
```

rewrites all five figures into `figures/` from `data/` alone.

## Software dependencies

| Domain | Package | Used by |
| --- | --- | --- |
| Residue code map and BLOSUM62 | Biopython | 2 of the deposited scripts |
| Figures | Matplotlib | 1 of the deposited scripts |
| Arrays and statistics | NumPy | 11 of the deposited scripts |
| Optimisation and spatial routines | SciPy | 1 of the deposited scripts |

Everything else the scripts use is in the Python standard library.

## Data availability

Four primary sources are used and none is redistributed here:

* cDNA-proteolysis folding stabilities (Tsuboyama et al. 2023), the FireProtDB curated set and the S669 benchmark, read from the distributions released with the ThermoMPNN repository.
* AlphaMissense scores over a 66-protein soluble cohort and a 504-sequence transmembrane cohort, from the published release.
* UniProt humsavar, release 2026_03.
* Deposited structures for the cDNA-proteolysis identifiers, obtained from the RCSB PDB by `codes/5-geometry/fetch_structures.py`, which reads `data/ladder_geometry_structures.csv` and checks each file on arrival.

## Citation

Karagöl, T., & Karagöl, A. (2026). Directional cost of an engineered histidine in near-isosteric charge-altering substitutions.

## License

The code in this repository is released under the MIT License. The upstream data follow the licence of each provider named above.
