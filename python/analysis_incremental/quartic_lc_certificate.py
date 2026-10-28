"""Exact temporal quartic certificate for the m=3 series-LC example.

Theorem 3(a) of the manuscript (Theorem 3 of the extended version).
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



# Purely temporal input: u=sqrt(2)*P(cos(t)), i=G*u+h(A)*Du.
P = [-1, -2, 2, 4]
p = [-36, 49, -14, 1]
h = [Fraction(-322,180), Fraction(155,180), Fraction(-13,180)]
krylov = [deriv(P)]
for _ in range(3):
    krylov.append(temporal_A(krylov[-1]))
assert reduce(add, [scale(q, pj) for q,pj in zip(krylov,p)], [0]) == [0]
K = [[krylov[j][k] if k<len(krylov[j]) else 0 for j in range(3)] for k in range(3)]
rank_det = determinant(K)
assert rank_det != 0
it = reduce(add, [scale(q,hj) for q,hj in zip(krylov,h)], [0])
W = [[0] for _ in range(3)]
for j in range(1,len(p)):
    for k in range(j):
        W[j-1-k] = add(W[j-1-k], scale(krylov[k],p[j]))
a = [temporal_dot(it,q) for q in W]
b = [Fraction(0)]*5
for j,q in enumerate(W):
    for k,r in enumerate(W):
        b[j+k] += temporal_dot(q,r)
a,b = [int(x) for x in a], [int(x) for x in b]
assert a == [-14,-12,2]
assert b == [1764,-2016,1010,-196,14]
raw = add(scale(mul(deriv(a),b),2), scale(mul(a,deriv(b)),-1))
H = scale(primitive(raw),-1)
assert H == [8820,-8323,2541,-309,7]
assert raw == scale(H,-8)
# The prescribed load current is realized by a fully passive parallel RLC load.
# h(A) = 11/20 - (17/10)/A + (3/4)/(6-A), modulo p(A)=0.
# A(6-A) is invertible on the cyclic space because p(0),p(6) are nonzero.
A_six_minus_A=[0,6,-1]
passive_difference=add(
    add(mul(h,A_six_minus_A), scale(A_six_minus_A,Fraction(-11,20))),
    add(scale([6,-1],Fraction(17,10)), scale([0,1],Fraction(-3,4))))
passive_identity_rhs=scale(mul([-51,13],p),Fraction(1,180))
assert passive_difference==passive_identity_rhs
assert evaluate(p,0)!=0 and evaluate(p,6)!=0
passive_components={'shunt_capacitance':Fraction(11,20),
                    'shunt_inductance':Fraction(10,17),
                    'series_branch_inductance':Fraction(4,3),
                    'series_branch_capacitance':Fraction(1,8)}
assert all(value>0 for value in passive_components.values())
assert 1/passive_components['series_branch_inductance']==Fraction(3,4)
assert 1/(passive_components['series_branch_inductance']*
           passive_components['series_branch_capacitance'])==6
# Rabin test for irreducibility of a quartic: gcd(f,x^(p^2)-x)=1
# and x^(p^4)=x modulo f.
prime=13
Hp=modp(H,prime)
x=[0,1]
r=x
for _ in range(2):
    r=powmodp(r,prime,Hp,prime)
g=gcdmodp(Hp,modp(add(r,scale(x,-1)),prime),prime)
r2=r
for _ in range(2):
    r=powmodp(r,prime,Hp,prime)
assert g == [1]
assert r == x
# Optional strengthening: a 4-cycle modulo 13 and a 3-cycle modulo 19
# certify S4, so square roots alone also do not suffice, while general
# radicals (Ferrari) still do.
linear19=[18,1]
cubic19=[15,16,2,7]
assert modp(mul(linear19,cubic19),19)==modp(H,19)
assert irreducible_prime_degree(cubic19,19)
assert len(gcdmodp(modp(H,19),deriv(modp(H,19)),19))==1
poles=[1,4,9]
pole_gains=[Fraction(evaluate(a,z)**2,evaluate(b,z)) for z in poles]
gain2=Fraction(evaluate(a,2)**2,evaluate(b,2))
assert gain2 == Fraction(225,107)
assert gain2 > max(pole_gains)
assert evaluate(p,0)*evaluate(a,0)>0  # inactive L-only boundary
assert a[-1]>0                       # inactive C-only boundary
chain=sturm(H)
bound=Fraction(1+max(abs(v) for v in H[:-1])//abs(H[-1])+1)
segments=[(Fraction(lo),Fraction(hi)) for lo,hi in zip([0,1,4,9],[1,4,9,bound])]
isolated=[]
while segments:
    lo,hi=segments.pop()
    count=variations(chain,lo)-variations(chain,hi)
    if not count:
        continue
    if count==1:
        isolated.append((lo,hi))
        continue
    mid=(lo+hi)/2
    segments.extend([(lo,mid),(mid,hi)])
roots=[]
for lo,hi in sorted(isolated):
    for _ in range(80):
        mid=(lo+hi)/2
        if variations(chain,lo)-variations(chain,mid):
            hi=mid
        else:
            lo=mid
    z=float((lo+hi)/2)
    av,bv,pv=evaluate(a,z),evaluate(b,z),evaluate(p,z)
    active=pv*av<0
    gamma=max(0.,-pv*av/bv)
    roots.append({'z':z,'active':active,'gain':av*av/bv if active else 0.,
                  'gamma':gamma,'L':1/gamma if gamma else None,
                  'C':gamma/z if gamma else None})
Gram=[[temporal_dot(q,r) for r in krylov[:3]] for q in krylov[:3]]
correlations=[temporal_dot(it,q) for q in krylov[:3]]
assert Gram==[[14,98,794],[98,794,6818],[794,6818,60074]]
assert correlations==[2,16,112]
stationary_counts=[]
for lo,hi in zip([0,1,4,7,9],[1,4,7,9,bound]):
    count=variations(chain,Fraction(lo))-variations(chain,Fraction(hi))
    stationary_counts.append({'lo':str(lo),'hi':str(hi),'count':count})
assert [r['count'] for r in stationary_counts]==[0,1,0,0,1]
isolating_lo,isolating_hi=Fraction(1986,1000),Fraction(1987,1000)
assert variations(chain,isolating_lo)-variations(chain,isolating_hi)==1
assert len([r for r in roots if r['active']])==1
assert all(r['L']>0 and r['C']>0 for r in roots if r['active'])
result={
    'time_input': {'u':'sqrt(2)*(4*cos(t)^3+2*cos(t)^2-2*cos(t)-1)',
                   'i':'G*u+h(A)*Du; A=-D^2; G arbitrary nonnegative',
                   'h_ascending_exact':[str(x) for x in h]},
    'explicit_passive_load_certificate': {
        'topology':'Parallel conductance G, capacitor C0, inductor L0, and one series LC branch',
        'components_exact':{key:str(value) for key,value in passive_components.items()},
        'conductance':'Any G>=0; G=1 gives an explicitly dissipative passive load',
        'current_equation':'i=G*u+(11/20)*D*u+(17/10)*D^(-1)*u+(3/4)*(D^2+6)^(-1)*D*u',
        'rational_operator_identity':'h(A)=11/20-(17/10)*A^(-1)+(3/4)*(6-A)^(-1) modulo p(A)',
        'cleared_denominator_identity':'A*(6-A)*(h(A)-11/20+(17/10)*A^(-1)-(3/4)*(6-A)^(-1))=(13*A-51)*p(A)/180',
        'cleared_difference_ascending_exact':[str(value) for value in passive_difference],
        'identity_verified_exactly':True,
        'denominator_invertibility':{'p_at_0':evaluate(p,0),'p_at_6':evaluate(p,6)},
        'all_components_strictly_positive':True,
        'meaning':'The irrational optimal compensator tuning occurs even for an explicit passive load, independently of any universal synthesis theorem.'},
    'temporal_certificate': {'krylov_polynomials_ascending':krylov,
                             'rank_determinant_exact':str(rank_det),
                             'p_annihilation_verified':True,
                             'current_polynomial_ascending_exact':[str(x) for x in it],
                             'p_ascending':p,'a_ascending':a,'b_ascending':b,
                             'H_ascending':H,'F_equals':'-8*H',
                             'GA_basis':['Du','A*Du','A^2*Du'],
                             'GA_basis_Gram_exact':[[str(v) for v in row] for row in Gram],
                             'current_correlations_exact':[str(v) for v in correlations],
                             'u_dot_each_GA_basis_vector':[0,0,0],
                             'u_norm_squared_exact':'3',
                             'i_Q_norm_squared_exact':str(temporal_dot(it,it))},
    'irreducibility_certificate': {'prime':13,'H_mod_13_ascending':Hp,
                                   'x_to_p_squared_remainder':r2,
                                   'gcd_with_x_to_p_squared_minus_x':g,
                                   'x_to_p_fourth_remainder':r,
                                   'conclusion':'H irreducible over Q; each optimal tuning has algebraic degree 4'},
    'optional_galois_strengthening':{'group':'S4',
                                     'mod_13_cycle_type':[4],
                                     'mod_19_cycle_type':[1,3],
                                     'mod_19_factors_ascending':[linear19,cubic19],
                                     'cubic_irreducible_verified':True,
                                     'meaning':'Products, quotients and nested square roots alone do not suffice; general radicals still suffice for this quartic.'},
    'interior_existence': {'gain_at_z_2_exact':str(gain2),
                           'gain_pole_limits_exact':[str(x) for x in pole_gains],
                           'gain_z_0':0,'gain_z_infinity':0},
    'positive_stationary_roots':roots,
    'uniqueness_certificate':{'active_intervals':['(1,4)','(7,9)'],
                               'stationary_root_counts':stationary_counts,
                               'optimum_isolating_interval_exact':[str(isolating_lo),str(isolating_hi)],
                               'unique_active_stationary_root':True,
                               'reason':'The sole active stationary root is in (1,4), where z>0 and p*a<0 guarantee Gamma,L,C>0. Its gain exceeds all boundary limits. The other active interval has no stationary root. Hence the positive finite LC minimizer is globally unique.'},
    'global_optimum':max(roots,key=lambda r:r['gain'])
}
destination=Path(__file__).with_suffix('.json')
destination.write_text(json.dumps(result,indent=2)+'\n')
print('Quartic:',H)
print('Irreducible modulo 13, independently verified by Rabin test.')
print('Interior gain:',gain2,'; all excluded-boundary gains at most',max(pole_gains))
print('Positive stationary roots:',roots)
print('Global optimum:',result['global_optimum'])
print('Explicit passive RLC load identity verified exactly; all component values positive.')
print('Certificate saved:',destination)
