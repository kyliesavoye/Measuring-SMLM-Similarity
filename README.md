# Measuring-SMLM-Similarity

Code to generate the simulations, analyses, and figures for:

> Savoye, K., Nieves, D.J., Shirgill, S., Lewis, A., Spill, F., Owen, D.M. (2025).
> **Measuring the similarity of single-molecule localization microscopy derived marked point clouds.**
> *Biophysical Journal*, 124, 2931-2940. https://doi.org/10.1016/j.bpj.2025.07.035

This repository is a maintained, updated re-publication of the original project.

## Overview

Single-molecule localization microscopy (SMLM) combined with environmentally sensitive fluorescent probes (e.g. di-4-ANEPPDHQ) produces **marked point clouds**: sets of spatial coordinates each carrying an additional value (a "mark"), such as a generalized polarization (GP) value reporting on membrane lipid order.

This repository implements a method that extends the two-sample Kolmogorov-Smirnov (KS) test to marked point clouds. For any pair of point clouds, it computes three semi-independent KS scores:

1. **Coordinate KS score** - compares the spatial distribution of localizations (grid-binned cumulative histograms), ignoring marks.
2. **Mark (GP) KS score** - compares the cumulative distribution of mark values directly, ignoring spatial position.
3. **Sum-of-marks KS score** - compares cumulative histograms of the *summed* mark values per spatial bin, capturing how marks are spatially arranged.

The three scores define a point in 3D space for each pairwise comparison; the Euclidean distance to the origin, **L**, is used as a single similarity metric (smaller L = more similar).

The method is validated on simulated marked point clouds with varied spatial and mark-distribution parameters, then applied to experimental SMLM data comparing control vs. methyl-beta-cyclodextrin (MBCD)-treated cells, where it detects condition-dependent changes in membrane lipid order.

## Repository Contents

| Path | Description |
|---|---|
| `smlm/` | Core package: functions for computing the three KS scores, simulating marked point cloud data, and calculating the distance-to-origin similarity metric. |
| `Measuring_the_similarity_of_SMLM_marked_point_clouds.ipynb` | Main notebook reproducing the analysis pipeline end-to-end, from simulation through to the experimental comparison. |
| `Method Schematic Figues/` | Figures for the method schematic (Fig. 1). |
| `Simulation Data Example Figures/` | Figures for simulated data examples and parameter-variation figures (Figs. 2-5). |
| `Testing Method Figures/` | Figures for the parameter-testing and validation figures. |
| `Experimental Data/` | Experimental SMLM/GP data (control vs. MBCD-treated ROIs). |
| `Experimental Data Figure/` | Figures for the experimental comparison figure (Fig. 6). |

## Method Summary

For two samples drawn from distributions F and G, the two-sample KS statistic is the maximum distance between their empirical cumulative distribution functions. This is rescaled into a KS score (Eq. 1 in the paper) such that a score of 0 indicates identical distributions, with higher scores indicating greater dissimilarity.

The three KS scores (coordinate, mark, sum-of-marks) are computed per comparison and combined into a single distance-to-origin metric:

```text
L = sqrt(KS1^2 + KS2^2 + KS3^2)
```

See the paper's Materials and Methods for full derivations, and the "Testing the method" and "Application" sections for how the metric is validated on simulated data and applied to experimental GP-marked SMLM data.

## Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/kyliesavoye/Measuring-SMLM-Similarity.git
cd Measuring-SMLM-Similarity
pip install -r requirements.txt
```

## Usage

Open `Measuring_the_similarity_of_SMLM_marked_point_clouds.ipynb` to reproduce the simulation studies and experimental comparison.

## Data

Experimental data (control and MBCD-treated ROIs) originate from the SMLM dataset described in Panconi et al., *Nat. Commun.* 15:9641 (2024), and are provided in `Experimental Data/` for reproducing the paper's figures.

Source dataset citation:

```bibtex
@article{panconi2024mapping,
	title={Mapping membrane biophysical nano-environments},
	author={Panconi, Luca and Euchner, Jonas and Tashev, Stanimir A and Makarova, Maria and Herten, Dirk-Peter and Owen, Dylan M and Nieves, Daniel J},
	journal={Nature Communications},
	volume={15},
	number={1},
	pages={9641},
	year={2024},
	publisher={Nature Publishing Group UK London}
}
```

## Citation

If you use this code, please cite:

```bibtex
@article{savoye2025measuring,
	title   = {Measuring the similarity of single-molecule localization microscopy derived marked point clouds},
	author  = {Savoye, Kylie and Nieves, Daniel J. and Shirgill, Sandeep and Lewis, Arthur and Spill, Fabian and Owen, Dylan M.},
	journal = {Biophysical Journal},
	volume  = {124},
	pages   = {2931--2940},
	year    = {2025},
	doi     = {10.1016/j.bpj.2025.07.035}
}
```

## Funding

This work was supported by the EPSRC Centre for Doctoral Training in Topological Design (EP/S02297X/1, K.S.), AstraZeneca, a UKRI Future Leaders Fellowship (MR/T043571/1, F.S.), and BBSRC (BB/X018644/1).

## License

See [LICENSE](LICENSE).
