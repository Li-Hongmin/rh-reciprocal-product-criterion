# Reproducing the finite checks

These files accompany *A Positive Reciprocal-Product Criterion for the
Riemann Hypothesis*, submitted by the account `lihongmin`.

Run from this directory:

```sh
python3 verify_reciprocal_kernel.py
uv run --with sympy python verify_detector_certificates.py
```

The first command uses the Python standard library and checks finite kernel
identities and two squarefree-star examples in exact rational arithmetic.
The second reconstructs the five all-half-line detector certificates and
the five minimality ranks with SymPy. Expected output: two PASS lines from
the first command; five detector PASS lines and one rank PASS line from
the second. The manuscript supplies the analytic arguments and their scope.

The file supplementary-detector-certificates.md contains the exact tables.
No numerical sample is used as evidence of RH or an infinite hierarchy.
The verification workload is a single CPU process, not a large simulation.
See LICENSE.txt for rights and permitted verification use.

- `verify_degree_bound.py`: equality case (8.0b) of Proposition 8.1 for the five detectors; index condition of Corollary 6.1. Run `uv run --with sympy python verify_degree_bound.py`.
