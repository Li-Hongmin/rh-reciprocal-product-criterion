"""Exact symbolic checks for the finite detector results in the manuscript.

Run with:
    uv run --with sympy python verify_detector_certificates.py
"""

from functools import lru_cache

import sympy as sp


u, s, N, r, m, t, c = sp.symbols("u s N r m t c")


DETECTORS = {
    1: u**2 * (1 - u) ** 3 * (7 * u - 3),
    2: u**2 * (1 - u) ** 5 * (20 * u**2 - 5 * u - 1),
    3: u**2 * (1 - u) ** 7 * (91 * u**3 - 3 * u**2 - 7 * u - 1),
    4: u**2
    * (1 - u) ** 9
    * (1024 * u**4 + 207 * u**3 - 57 * u**2 - 27 * u - 3),
    5: u**2
    * (1 - u) ** 11
    * (8151 * u**5 + 3645 * u**4 + 198 * u**3 - 262 * u**2 - 77 * u - 7)
    / 8151,
}


EXPECTED_MELLIN = {
    1: 24 * (s - 1) / sp.prod(s + k for k in (2, 3, 4, 5, 6)),
    2: 1680 * (s - 1) / sp.prod(s + k for k in (2, 4, 5, 6, 7, 8, 9)),
    3: 403200
    * (s - 1)
    / sp.prod(s + k for k in (2, 4, 6, 7, 8, 9, 10, 11, 12)),
    4: 415134720
    * (s - 1)
    / sp.prod(s + k for k in (2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15)),
    5: 1083801600
    * (s - 1)
    / (
        19
        * sp.prod(
            s + k for k in (2, 4, 6, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18)
        )
    ),
}


@lru_cache(maxsize=None)
def power_sum(k):
    return sp.expand(sp.summation(c**k, (c, 1, N)))


def mellin_transform(poly):
    result = 0
    for (degree,), coefficient in sp.Poly(sp.expand(poly), u).terms():
        result += coefficient / (s + degree)
    return sp.factor(result)


def bernstein_coefficients(poly):
    """Return coefficients of Q_N in the degree-(d-2) Bernstein basis.

    P_V(N+r) = N(N+1) Q_N(r) / (N+r)^d.
    """

    degree = sp.degree(poly, u)
    summed = 0
    for (k,), coefficient in sp.Poly(sp.expand(poly), u).terms():
        summed += coefficient * power_sum(k) / (N + r) ** k
    q_poly = sp.Poly(
        sp.expand(sp.cancel(summed * (N + r) ** degree / (N * (N + 1)))), r
    )
    bernstein_degree = degree - 2
    powers = [
        q_poly.coeff_monomial(r**j) for j in range(bernstein_degree + 1)
    ]
    return [
        sp.factor(
            sum(
                powers[j]
                * sp.binomial(k, j)
                / sp.binomial(bernstein_degree, j)
                for j in range(k + 1)
            )
        )
        for k in range(bernstein_degree + 1)
    ]


def falling_coefficients(poly, variable, degree):
    """Coefficients in 1, x, x(x-1), ... by finite differences."""

    expanded = sp.expand(poly)
    return [
        sp.factor(
            sum(
                (-1) ** (j - i) * sp.binomial(j, i) * expanded.subs(variable, i)
                for i in range(j + 1)
            )
            / sp.factorial(j)
        )
        for j in range(degree + 1)
    ]


def assert_nonnegative_rationals(values, label):
    for value in values:
        value = sp.factor(value)
        assert value.is_Rational and value >= 0, (label, value)


def verify_detector(level, poly):
    poly = sp.expand(poly)
    degree = 3 * level + 3
    assert sp.degree(poly, u) == degree
    assert sp.simplify(mellin_transform(poly) - EXPECTED_MELLIN[level]) == 0
    assert sp.integrate(poly, (u, 0, 1)) == 0
    for odd_degree in range(1, 2 * level, 2):
        assert sp.Poly(poly, u).coeff_monomial(u**odd_degree) == 0

    bernstein = bernstein_coefficients(poly)
    if level <= 4:
        for k, coefficient in enumerate(bernstein):
            values = falling_coefficients(
                coefficient.subs(N, m + 1), m, level
            )
            assert_nonnegative_rationals(values, f"V{degree}, k={k}")
    else:
        at_one = [sp.factor(value.subs(N, 1)) for value in bernstein]
        assert_nonnegative_rationals(at_one, "V18, N=1")
        for k, coefficient in enumerate(bernstein[1:], start=1):
            polynomial = sp.Poly(sp.expand(coefficient.subs(N, t + 2)), t)
            assert_nonnegative_rationals(
                polynomial.all_coeffs(), f"V18, N>=2, k={k}"
            )
        expected_b0 = (
            (N - 1)
            * (1746848 * N**4 - 14791632 * N**2 + 31364905)
            / 570570
        )
        assert sp.simplify(bernstein[0] - expected_b0) == 0
        special_factor = 1746848 * N**4 - 14791632 * N**2 + 31364905
        assert sp.expand(
            54589 * special_factor
            - 2 * (218356 * N**2 - 924477) ** 2
        ) == 2863351987


def verify_minimality_ranks():
    for level in range(1, 6):
        coefficients = sp.symbols(f"a0:{level}")
        quotient = sum(coefficients[j] * u**j for j in range(level))
        candidate = sp.expand(u**2 * (1 - u) ** (2 * level + 1) * quotient)
        equations = [
            sp.Poly(candidate, u).coeff_monomial(u ** (2 * j + 1))
            for j in range(1, level)
        ]
        equations.append(sp.integrate(candidate, (u, 0, 1)))
        matrix, _ = sp.linear_eq_to_matrix(equations, coefficients)
        assert matrix.rank() == level


def main():
    for level, detector in DETECTORS.items():
        verify_detector(level, detector)
        print(f"V_{3 * level + 3}: PASS")
    verify_minimality_ranks()
    print("minimality ranks m=1,...,5: PASS")


if __name__ == "__main__":
    main()
