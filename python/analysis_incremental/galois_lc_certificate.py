"""Exact finite-field and real-root certificate for the m=4 series-LC example.

Theorem 3(b) of the manuscript (Theorem 4 of the extended version).
No third-party packages. Polynomials use ascending coefficients.
"""
from fractions import Fraction
from functools import reduce
from math import gcd, lcm, comb
from pathlib import Path
import json


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def add(a, b):
    return trim([(a[k] if k < len(a) else 0) +
                 (b[k] if k < len(b) else 0)
                 for k in range(max(len(a), len(b)))])


def scale(a, t):
    return trim([x*t for x in a])


def mul(a, b):
    c = [0]*(len(a)+len(b)-1)
    for j, x in enumerate(a):
        for k, y in enumerate(b):
            c[j+k] += x*y
    return trim(c)


def deriv(a):
    return [j*a[j] for j in range(1, len(a))] or [0]


def evaluate(a, x):
    y = 0
    for c in reversed(a):
        y = y*x+c
    return y


def modp(a, prime):
    return trim([x % prime for x in a])


def divmodp(a, b, prime):
    a, b = modp(a, prime), modp(b, prime)
    q = [0]*max(1, len(a)-len(b)+1)
    while a != [0] and len(a) >= len(b):
        k = len(a)-len(b)
        c = a[-1]*pow(b[-1], -1, prime) % prime
        q[k] = c
        for j, v in enumerate(b):
            a[k+j] = (a[k+j]-c*v) % prime
        a = trim(a)
    return trim(q), a


def gcdmodp(a, b, prime):
    while b != [0]:
        a, b = b, divmodp(a, b, prime)[1]
    return modp(scale(a, pow(a[-1], -1, prime)), prime)


def powmodp(a, n, f, prime):
    result = [1]
    while n:
        if n & 1:
            result = divmodp(mul(result, a), f, prime)[1]
        n //= 2
        if n:
            a = divmodp(mul(a, a), f, prime)[1]
    return result


def irreducible_prime_degree(f, prime):
    """Rabin test specialized to prime degrees 2, 3, 5, 7."""
    degree = len(f)-1
    assert degree in (2, 3, 5, 7)
    x = [0, 1]
    h = x
    for _ in range(degree):
        h = powmodp(h, prime, f, prime)
    if modp(add(h, scale(x, -1)), prime) != [0]:
        return False
    g = gcdmodp(f, modp(add(powmodp(x, prime, f, prime),
                              scale(x, -1)), prime), prime)
    return len(g) == 1


def primitive(a):
    a = trim([Fraction(v) for v in a])
    den = reduce(lcm, (v.denominator for v in a), 1)
    ints = [int(v*den) for v in a]
    content = reduce(gcd, ints)
    return [v//content for v in ints] if content else [0]


def rational_remainder(a, b):
    a, b = [Fraction(v) for v in a], [Fraction(v) for v in b]
    while a != [0] and len(a) >= len(b):
        k = len(a)-len(b)
        c = a[-1]/b[-1]
        for j, v in enumerate(b):
            a[k+j] -= c*v
        a = trim(a)
    return a


def sturm(a):
    chain = [primitive(a), primitive(deriv(a))]
    while len(chain[-1]) > 1:
        nxt = primitive(scale(rational_remainder(chain[-2], chain[-1]), -1))
        assert nxt != [0]
        chain.append(nxt)
    return chain


def variations(chain, x):
    values = [evaluate(a, x) for a in chain]
    signs = [1 if v > 0 else -1 for v in values if v]
    return sum(a != b for a, b in zip(signs, signs[1:]))


def temporal_A(q):
    """A on -sqrt(2)*sin(t)*q(cos(t)), obtained by two differentiations."""
    return add(add(q, scale(mul([0, 1], deriv(q)), 3)),
               mul([-1, 0, 1], deriv(deriv(q))))


def cos_moment(power):
    if power % 2:
        return Fraction(0)
    return Fraction(comb(power, power//2), 2**power)


def temporal_dot(q, r):
    """Exact period integral of 2*sin(t)^2*q(cos(t))*r(cos(t))."""
    return sum(Fraction(2)*coef*(cos_moment(k)-cos_moment(k+2))
               for k, coef in enumerate(mul(q, r)))


def determinant(matrix):
    matrix = [[Fraction(x) for x in row] for row in matrix]
    result = Fraction(1)
    for j in range(len(matrix)):
        pivot = next((k for k in range(j, len(matrix)) if matrix[k][j]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != j:
            matrix[j], matrix[pivot] = matrix[pivot], matrix[j]
            result = -result
        value = matrix[j][j]
        result *= value
        for k in range(j+1, len(matrix)):
            ratio = matrix[k][j]/value
            matrix[k] = [x-ratio*y for x, y in zip(matrix[k], matrix[j])]
    return result


lam = [1, 4, 9, 16]
d = [1, 4, 9, 16]
c = [-7, 14, 9, 9]
p = reduce(mul, [[-x, 1] for x in lam], [1])
quotients = [reduce(mul, [[-x, 1] for j, x in enumerate(lam) if j != k], [1])
             for k in range(4)]
a = reduce(add, [scale(q, cj) for q, cj in zip(quotients, c)], [0])
b = reduce(add, [scale(mul(q, q), dj) for q, dj in zip(quotients, d)], [0])
F = primitive(add(scale(mul(deriv(a), b), 2),
                  scale(mul(a, deriv(b)), -1)))
h = reduce(add, [scale(q, Fraction(cj, dj*evaluate(deriv(p), lj)))
                 for q, cj, dj, lj in zip(quotients, c, d, lam)], [0])
assert all(evaluate(h, lj) == Fraction(cj, dj)
           for cj, dj, lj in zip(c, d, lam))

# Independent temporal construction: no decomposition into harmonic components.
# u(t)=sqrt(2)*P(cos(t)); Du=-sqrt(2)*sin(t)*P'(cos(t)).
P_time = [0, -2, -6, 4, 8]
krylov_time = [deriv(P_time)]
for _ in range(4):
    krylov_time.append(temporal_A(krylov_time[-1]))
assert reduce(add, [scale(q, pj) for q, pj in zip(krylov_time, p)], [0]) == [0]
coefficient_matrix = [[krylov_time[j][k] if k < len(krylov_time[j]) else 0
                       for j in range(4)] for k in range(4)]
krylov_determinant = determinant(coefficient_matrix)
assert krylov_determinant != 0
i_time = reduce(add, [scale(q, hj) for q, hj in zip(krylov_time, h)], [0])
W_time = [[0] for _ in range(4)]
for j in range(1, len(p)):
    for k in range(j):
        W_time[j-1-k] = add(W_time[j-1-k], scale(krylov_time[k], p[j]))
a_time = [temporal_dot(i_time, q) for q in W_time]
b_time = [Fraction(0) for _ in range(7)]
for j, q in enumerate(W_time):
    for k, r in enumerate(W_time):
        b_time[j+k] += temporal_dot(q, r)
assert a_time == a
assert b_time == b
assert temporal_dot(i_time, i_time) == Fraction(1793, 16)

assert F == [68173952, -80912600, 36926304, -8601960,
             1041441, -60975, 1003, 35]
assert irreducible_prime_degree(modp(F, 31), 31)
quadratic = [28, 14, 1]
quintic = [2, 18, 27, 28, 20, 6]
assert modp(mul(quadratic, quintic), 29) == modp(F, 29)
assert irreducible_prime_degree(quadratic, 29)
assert irreducible_prime_degree(quintic, 29)
assert len(gcdmodp(modp(F, 29), deriv(modp(F, 29)), 29)) == 1

print('Exact S7 certificates verified: irreducible mod 31; degrees 2+5 mod 29.')
print('Ascending coefficients:', {'p': p, 'a': a, 'b': b, 'F': F})
print('Temporal current i_Q=h(A)Du; h coefficients ascending:', [str(v) for v in h])

gain_at_two = evaluate(a, 2)**2 / Fraction(evaluate(b, 2))
pole_gains = [Fraction(cj*cj, dj) for cj, dj in zip(c, d)]
assert gain_at_two == Fraction(49729, 444)
assert gain_at_two > max(pole_gains)
assert sum(c) > 0  # C-only inactive.
assert sum(Fraction(cj, lj) for cj, lj in zip(c, lam)) < 0  # L-only inactive.
print('Interior feasibility:', {'gain_at_2': str(gain_at_two),
                               'pole_gain_max': str(max(pole_gains)),
                               'L_C_boundary_gains': 0})

chain = sturm(F)
bound = Fraction(1 + max(abs(v) for v in F[:-1]) // abs(F[-1]) + 1)
segments = [(Fraction(lo), Fraction(hi)) for lo, hi in
            zip([0, 1, 4, 9, 16], [1, 4, 9, 16, bound])]
isolated = []
while segments:
    lo, hi = segments.pop()
    count = variations(chain, lo)-variations(chain, hi)
    if count == 0:
        continue
    if count == 1:
        isolated.append((lo, hi))
        continue
    middle = (lo+hi)/2
    assert evaluate(F, middle) != 0  # Irreducibility excludes rational roots.
    segments.extend([(lo, middle), (middle, hi)])

roots = []
for lo, hi in sorted(isolated):
    for _ in range(80):
        mid = (lo+hi)/2
        if variations(chain, lo)-variations(chain, mid):
            hi = mid
        else:
            lo = mid
    z = float((lo+hi)/2)
    av, bv, pv = evaluate(a, z), evaluate(b, z), evaluate(p, z)
    active = pv*av < 0
    gamma = max(0., -pv*av/bv)
    gain = av*av/bv if active else 0.
    roots.append({'z': z, 'active': active, 'gain': gain, 'gamma': gamma,
                  'L': 1/gamma if gamma else None,
                  'C': gamma/z if gamma else None})
print('All positive stationary roots isolated by exact Sturm counts:', roots)
print('Global optimum:', max(roots, key=lambda r: r['gain']))


def rabin_details(poly, prime):
    x = [0, 1]
    degree = len(poly)-1
    h = x
    for _ in range(degree):
        h = powmodp(h, prime, poly, prime)
    first_power_difference = modp(add(powmodp(x, prime, poly, prime),
                                     scale(x, -1)), prime)
    return {'prime': prime, 'degree': degree,
            'polynomial_ascending': modp(poly, prime),
            'x_to_p_to_degree_remainder': h,
            'gcd_with_x_to_p_minus_x': gcdmodp(poly, first_power_difference, prime),
            'irreducible': irreducible_prime_degree(poly, prime)}


certificate = {
    'formulation': {
        'u_time': 'sqrt(2)*(8*cos(t)**4+4*cos(t)**3-6*cos(t)**2-2*cos(t))',
        'T': '2*pi', 'A': '-D**2', 'v0': 'D*u',
        'i': 'G*u+h(A)*D*u, with arbitrary G>=0',
        'h_ascending_exact': [str(v) for v in h],
        'p_ascending': p, 'a_ascending': a, 'b_ascending': b,
        'F_primitive_ascending': F,
        'stationary_polynomial_scale': '2*a_prime*b-a*b_prime = 36*F',
        'auxiliary_eigenvalues_for_certificate_only': lam,
        'auxiliary_d': d, 'auxiliary_c': c,
    },
    'finite_field_certificates': {
        'irreducible_mod_31': rabin_details(F, 31),
        'factorization_mod_29': {
            'factors_ascending': [quadratic, quintic],
            'product_verified': modp(mul(quadratic, quintic), 29) == modp(F, 29),
            'squarefree_gcd': gcdmodp(modp(F, 29), deriv(modp(F, 29)), 29),
            'quadratic_rabin': rabin_details(quadratic, 29),
            'quintic_rabin': rabin_details(quintic, 29),
        },
        'galois_group': 'S7',
        'reason': 'Irreducibility gives transitivity of prime degree 7. The unramified '
                  'factorization 2+5 yields an element whose fifth power is a '
                  'transposition. A transitive prime-degree subgroup with a '
                  'transposition is S7.',
    },
    'independent_temporal_certificate': {
        'coordinate': 'x=cos(t)',
        'u_polynomial_ascending': P_time,
        'A_on_q': 'q+3*x*q_prime+(x*x-1)*q_double_prime',
        'waveform_represented_by_q': '-sqrt(2)*sin(t)*q(cos(t))',
        'krylov_vectors_polynomial_ascending': krylov_time,
        'krylov_rank_determinant_exact': str(krylov_determinant),
        'p_annihilation_verified': True,
        'current_polynomial_ascending_exact': [str(v) for v in i_time],
        'W_coefficient_vectors_polynomial_ascending': W_time,
        'a_from_exact_period_integrals': [str(v) for v in a_time],
        'b_from_exact_period_integrals': [str(v) for v in b_time],
        'i_quadrature_norm_squared_exact': '1793/16',
    },
    'global_existence_certificate': {
        'gain_at_z_2_exact': str(gain_at_two),
        'gain_pole_limits_exact': [str(v) for v in pole_gains],
        'gain_z_0': 0, 'gain_z_infinity': 0,
        'argument': 'The feasible interior gain at z=2 exceeds every excluded '
                    'boundary limit, so a global optimum is attained in the '
                    'admissible interior and must be a root of F.',
    },
    'sturm_sequence_ascending': chain,
    'positive_roots': roots,
    'global_optimum': max(roots, key=lambda r: r['gain']),
}
destination = Path(__file__).with_suffix('.json')
destination.write_text(json.dumps(certificate, indent=2)+'\n')
print('Certificate written:', destination)
