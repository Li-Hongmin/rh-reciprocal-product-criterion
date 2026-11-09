"""Small exact checks of the kernel identities; not a proof of RH."""
from fractions import Fraction as F
from math import isqrt, prod


def mobius(n):
    sign = 1
    p = 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            sign = -sign
            if n % p == 0:
                return 0
        p += 1
    return -sign if n > 1 else sign


def v(u):
    return -2 * u + 6 * u**2 - 4 * u**3


def kernel(y):
    n = int(y)
    r = y - n
    return n * (n + 1) * r * (1 - r) / y**3


def inertia(matrix):
    a = [row[:] for row in matrix]
    signs = []
    for k in range(len(a)):
        pivot = a[k][k]
        assert pivot != 0
        signs.append(1 if pivot > 0 else -1)
        for i in range(k + 1, len(a)):
            for j in range(i, len(a)):
                a[i][j] -= a[i][k] * a[j][k] / pivot
                a[j][i] = a[i][j]
    return signs.count(1), signs.count(-1), 0


def check():
    assert [mobius(n) for n in range(1, 13)] == [1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0]
    for x in map(F, ['1', '1.5', '2', '6.25', '9', '15.75', '30', '64', '100.25']):
        z = isqrt(int(x))
        for n in range(1, int(x) + 1):
            convolution = sum(mobius(a) * mobius(b)
                for a in range(1, z + 1) for b in range(1, z + 1)
                if n % (a * b) == 0)
            assert mobius(n) == 2 * (mobius(n) if n <= z else 0) - convolution
        assert sum(v(F(c) / x) for c in range(1, int(x) + 1)) == kernel(x)
        u = sum(mobius(n) * v(F(n) / x) for n in range(1, int(x) + 1))
        a = sum(mobius(n) * v(F(n) / x) for n in range(1, z + 1))
        t = sum(mobius(i) * mobius(j) * kernel(x / (i * j))
            for i in range(1, z + 1) for j in range(1, z + 1))
        assert u == 2 * a - t
    assert kernel(F(3, 2)) == F(4, 27)
    assert kernel(F(9, 4)) == F(8, 81)
    assert -kernel(F(3, 2))**2 == F(-16, 729)
    primes = [2, 3, 5, 7, 11, 13, 17]
    for m in [6, 7]:
        x = prod(primes[:m])
        indices = [2 * p for p in primes[1:m]]
        k = [[kernel(F(x, a * b)) for b in indices] for a in indices]
        c = [[k[i][j] - F(a * b, 6 * x) for j, b in enumerate(indices)]
            for i, a in enumerate(indices)]
        assert inertia(k) == inertia(c) == (1, m - 2, 0)
    print('exact truncation, Faulhaber identity, and transfer: PASS')
    print('indefinite principal block and squarefree stars m=6,7: PASS')


if __name__ == '__main__':
    check()
