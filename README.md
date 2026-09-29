# PFAS-mixture-PBTK-Shapley

PBTK model and Shapley/Owen process-attribution code for PFAS mixture toxicokinetics in zebrafish.

This repository contains the deterministic physiologically based toxicokinetic (PBTK) model and counterfactual attribution workflow used to quantify process contributions to changes in PFAS internal dose under mixture exposure.

## Scope

The code implements:

- PBTK simulation of PFOS, PFHxS, PFOA, and PFBS in zebrafish
- single- and mixture-exposure toxicokinetic parameterization
- counterfactual simulations of toxicokinetic process changes
- first-level Shapley attribution of:
  - F: whole-blood unbound fraction
  - H: apparent hepatobiliary transfer
  - N: free-normalized apparent non-fecal clearance
- nested Owen decomposition of the unbound-fraction contribution into:
  - F1: blood-tissue exchange
  - F2: free-dependent non-fecal elimination
- calculation of modeled blood AUC reductions and process contributions

The attribution framework satisfies:

- F + H + N = total modeled AUC reduction
- F1 + F2 = F

## Repository Structure

```text
PFAS-mixture-PBTK-Shapley/
├── R/
│   └── PBTK_Shapley.R
├── README.md
├── requirements-R.txt
├── CITATION.cff
├── LICENSE
└── .gitignore
