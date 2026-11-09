# A Positive Reciprocal-Product Criterion for the Riemann Hypothesis

[![verify](https://github.com/Li-Hongmin/rh-reciprocal-product-criterion/actions/workflows/verify.yml/badge.svg)](https://github.com/Li-Hongmin/rh-reciprocal-product-criterion/actions/workflows/verify.yml) [![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.24054950-blue)](https://doi.org/10.5281/zenodo.24054950)

**Hongmin Li** (Institute of Science Tokyo; The University of Tokyo) · ORCID [0000-0003-0228-0600](https://orcid.org/0000-0003-0228-0600)

Main result of the manuscript: the Riemann hypothesis is equivalent to T(x) ≪ x^{1/2+ε} for an explicit sum T built from a positive reciprocal-product kernel P, and to the corresponding dyadic mean-square bound; the kernel's spectral structure and detector certificates are computed exactly.

- Manuscript on rStation: <https://rstation.ai/problems/riemann_hypothesis/papers/42/v3> (version 3)
- Archived on Zenodo: concept DOI [10.5281/zenodo.24054950](https://doi.org/10.5281/zenodo.24054950) (all versions); this repository mirrors version 3, DOI [10.5281/zenodo.24054974](https://doi.org/10.5281/zenodo.24054974)

## Contents

- `paper/paper.md`, `paper/rh-reciprocal-product-criterion.pdf` — the manuscript (version 3).
- `verification/` — the supplement of that version, unpacked unchanged; its own README describes each file.
- `.github/workflows/verify.yml` — reruns every check below on each push.

## What the programs check

| Check | Command (run inside `verification/`) |
|---|---|
| Finite kernel identities and squarefree-star examples | `python verify_reciprocal_kernel.py` |
| Five detector certificates and minimality ranks | `python verify_detector_certificates.py` |
| Degree bound (8.0b) and Corollary 6.1 index condition | `python verify_degree_bound.py` |

The programs check the finite and exact computations used in the manuscript; the mathematical arguments are those written in the paper.

## Reproducing

```sh
cd verification
pip install sympy==1.14.0
python verify_reciprocal_kernel.py
python verify_detector_certificates.py
python verify_degree_bound.py
```

Requirements: Python 3.9+. Run the commands from a copy if you want to keep the recorded outputs untouched.

## Cite

```bibtex
@misc{li2026rhreciprocalproductcrite,
  author    = {Li, Hongmin},
  title     = {A Positive Reciprocal-Product Criterion for the Riemann Hypothesis},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.24054950},
  url       = {https://doi.org/10.5281/zenodo.24054950}
}
```

## License

Manuscript and data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Code: MIT (see `LICENSE`). For this repository these licences, granted by the author, replace the earlier notice in `verification/LICENSE.txt`.
