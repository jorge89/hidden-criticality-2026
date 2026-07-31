Th# Hidden criticality

Code and data to reproduce the model results in:

* Fontenele, A.J. et al. Is critical brain dynamics more prevalent than previously thought? (submitted)

TL;DR: Simulates three network scenarios (Cases 1–3) and shows that near-critical dynamics can be hidden from the population-average signal, revealed only by analyzing individual modes.

## Usage

Run `compute_criticality_quantities.m`, setting `caseID` (1, 2, or 3) at the top, to simulate one of the three cases and compute per-mode criticality measures. Then run `makeFigs.m` to generate the corresponding panels for Fig. 2, Fig. 3, and Supplementary Fig. S2/S3.

d2 calculations depend on https://github.com/sam-sooter/prox_crit_toolkit.

(C) Antonio J. Fontenele and Woodrow L. Shew, 2026
