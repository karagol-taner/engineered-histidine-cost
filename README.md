# The cost of an engineered histidine

A prior on which substitutions a protein will accept when a histidine is installed, costed once against four readouts that share no measurements, and bounded by what it can and cannot be asked to decide.

This repository contains the source code and derived datasets for the analysis presented in the manuscript. It covers the per-substitution cost, the functional stratification, the within-protein and bootstrap robustness checks, the variance decomposition, the side-chain geometry, and every figure in the paper.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Citation

If you use this work in your research, please cite:

* Karagöl, T., & Karagöl, A. (2026). The cost of an engineered histidine: substitution tolerance replicates across measured folding stability, variant effect prediction and curated human variation.

## Usage

Every script takes `--help` and writes a CSV into `data/` plus a printed run log. The logs of the runs the manuscript quotes are in `data/source_runs/`.

* **`codes/1-cost/`**: `ladder_ddg.py` reads the three stability sets, the two AlphaMissense cohorts and the humsavar variant file, maps every stability column to a common sign, and writes the per-substitution Cliff delta for each readout.
* **`codes/2-stratification/`**: `ladder_humsavar_clusters.py` gives the curated-variation fractions a gene-clustered interval. `ladder_functional.py` splits a cohort into positions UniProt annotates as an active site, a binding site or a site and the rest, and scores each substitution against the baseline computed inside its own stratum.
* **`codes/3-robustness/`**: `ladder_robustness.py` ranks the substitutions inside each protein separately, resamples proteins rather than positions for the intervals, and stratifies by a leave-one-out constraint score. `ladder_designed.py` compares designed against natural domains.
* **`codes/4-variance/`**: `ladder_variance.py` separates the per-protein effect into a between-substitution and a protein-by-substitution component and predicts held-out domains. `ladder_extra.py` runs the descriptor regressions, the alternative baselines, the topology comparison and the transition counts.
* **`codes/5-geometry/`**: `ladder_overlap_all.py` widens the same measurement to every side-chain type and `ladder_overlap_stability.py` tests the overlap against measured folding cost over every substitution class. `ladder_geometry.py` measures the eight residue types in the deposited structures of the natural cDNA-proteolysis domains and writes the shape overlaps, reaches and representative side chains.
* **`codes/7-significance/`**: `ladder_significance.py` tests each substitution against its baseline on the cluster bootstrap and each rank agreement by permutation, with Benjamini-Hochberg applied across the twelve and across the eight comparisons.
* **`codes/6-figures/`**: `figures.py` rebuilds all five figures from `data/` alone, so they can be checked without repeating the analysis.
* **`data/`**: every derived table the manuscript draws a number from, plus the run logs. `cohort_soluble_66.csv` and `cohort_transmembrane_504.csv` list the UniProt accessions of the two AlphaMissense cohorts, so the scored sets are identified here and not only by reference. The soluble file also carries the `protein_index` used in `ladder_within_protein_soluble.csv` and `ladder_percell_soluble.csv`, both of which now name their proteins. `ladder_fireprot_positions.csv` holds the per-measurement folding cost, accessible surface area and functional flag that Supplementary Figure S1b is drawn from, so the figure rebuilds without the FireProtDB distribution.
* **`figures/`**: the five published figures at 600 dpi.
* **`Supplementary_Information.pdf`**: the supplementary tables and figure.

## Data

Four primary sources are used and none is redistributed here:

* cDNA-proteolysis folding stabilities (Tsuboyama et al. 2023), the FireProtDB curated set and the S669 benchmark, read from the distributions released with the ThermoMPNN repository.
* AlphaMissense scores over a 66-protein soluble cohort and a 504-sequence transmembrane cohort.
* UniProt humsavar, release 2026_03.
* Deposited structures for the 188 cDNA-proteolysis identifiers, downloaded from the RCSB PDB by the identifiers in `data/ladder_geometry_structures.csv`.

The substitution set and the single-nucleotide reachability test are in `codes/ladder_sets.py`. The set is the one the paper tabulates, and reachability is computed from the standard genetic code instead of being tabulated, so it can be read rather than trusted. The scripts use it by default and each one runs as deposited. `--targets` overrides it with `ladder_edge_cost.py` and `codon_audit.py`, which belong to the method of the citation above and are not distributed here; the two agree on every name the scripts use, and on reachability over all 400 ordered amino-acid pairs. The two published rankings are recomputed by running the method's own scripts when `--targets` is set, and read from their deposited raw output in `data/source_runs/` otherwise. `figures.py` uses none of this and rebuilds every figure from `data/` alone.

Each script shows its inputs under `--help`. The defaults expect the primary files beside the derived tables, as `data/megascale/`, `data/fireprot/`, `data/s669/`, `data/humsavar.txt`, `data/alphamissense_soluble.csv.gz` and `data/alphamissense_transmembrane/`, and `ladder_geometry.py` takes the structure directory explicitly. The derived tables in `data/` are the output of those runs, so every reported number and every figure can be checked without repeating them.

## Running the figures

```
cd codes/6-figures
python3 figures.py
```

The script reads only `data/` and writes into `figures/`. It locates the repository from its own path, or from `LADDER_ROOT` if that is set.

## Abstract

A protein that works in the bloodstream may have to switch off inside a tumour or an endosome, where the pH sits a unit lower. Engineers build that behaviour by installing histidine, the one amino acid that takes on charge across the physiological range, and they choose where to put it by screening, because nothing says in advance which positions will tolerate one. We costed every near-isosteric charge-altering substitution reachable by a single change at the first codon position, twelve in all, against four readouts that share no measurements: 271,231 experimental folding-stability measurements with two curated sets alongside, predicted pathogenicity over soluble and transmembrane cohorts of 306,622 and 5,846,493 substitutions, and 73,296 classified human variants. One ordering survives all four, at Spearman rho +0.622 to +0.909 across their six pairwise comparisons, and it reappears inside each of the 376 transmembrane proteins large enough to carry the test on its own. Converting an asparagine is the cheaper route to an installed histidine and converting an aspartate the dearer, by 35.8 against 55.4 per cent pathogenic variants and by +0.366 in Cliff delta. In a variance decomposition, which substitution it is carries 76 per cent of the effect a protein shows across all 361 substitution classes and 4 per cent of the variance at one position, so the ordering composes a scanning library and does not choose between two positions already characterised.

Keywords: histidine scanning, pH-responsive protein design, protein stability, substitution tolerance, variant effect prediction, pI-code
