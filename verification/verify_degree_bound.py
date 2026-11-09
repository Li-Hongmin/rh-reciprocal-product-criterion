"""Check that the five detectors V_6,...,V_18 attain the equality case (8.0b)
of Proposition 8.1: deg V = 3l+3 and K_V(s) = C(s-1)/(prod_{j<=l}(s+2j) prod_{k=2l+1}^{3l+3}(s+k)),
with C > 0 after positive normalization; also check the squarefree-index
condition 2p_m < sqrt(p_m#) used in Corollary 6.1 for m = 6..30."""
import sympy as sp
u, s = sp.symbols('u s')
V = {1: u**2*(1-u)**3*(7*u-3),
     2: u**2*(1-u)**5*(20*u**2-5*u-1),
     3: u**2*(1-u)**7*(91*u**3-3*u**2-7*u-1),
     4: u**2*(1-u)**9*(1024*u**4+207*u**3-57*u**2-27*u-3),
     5: u**2*(1-u)**11*(8151*u**5+3645*u**4+198*u**3-262*u**2-77*u-7)/8151}
for l, W in V.items():
    p = sp.Poly(sp.expand(W), u)
    assert p.degree() == 3*l+3
    assert all(p.coeff_monomial(u**(2*j-1)) == 0 for j in range(1, l+1))
    assert sp.integrate(W, (u, 0, 1)) == 0
    K = sp.together(sum(c/(s+m[0]) for m, c in zip(p.monoms(), p.coeffs())))
    D = sp.prod([s+2*j for j in range(1, l+1)]) * sp.prod([s+k for k in range(2*l+1, 3*l+4)])
    C = sp.simplify(K*D/(s-1))
    assert C.is_number and C > 0, (l, C)
    print(f"l={l}: deg {3*l+3}, K_V = {C}*(s-1)/D_l  OK")
P = 1
for m, q in enumerate(sp.primerange(2, 200), 1):
    P *= q
    if 6 <= m <= 30:
        assert 4*q*q < P
print("2p_m < sqrt(p_m#) for m=6..30 OK")
