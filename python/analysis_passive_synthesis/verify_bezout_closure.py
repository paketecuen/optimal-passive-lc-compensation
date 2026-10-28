"""Verification of temporal passive interpolation via Bezout matrices.

The synthesis routine receives p and h coefficients, not modal data.
Modal data below are ONLY used to generate/check test problems independently.
Requires NumPy and SciPy.
"""
from fractions import Fraction as F
import numpy as np
from numpy.polynomial import Polynomial as P
from scipy.linalg import eigh


def bezout(a, b):
    """B for (a(x)b(y)-a(y)b(x))/(x-y)=v(x)^T B v(y).

    a,b are numpy.polynomial.Polynomial (ascending coefficients).
    """
    n = max(a.degree(), b.degree())
    aa = np.pad(a.coef, (0, n + 1 - len(a.coef)))
    bb = np.pad(b.coef, (0, n + 1 - len(b.coef)))
    B = np.zeros((n, n))
    for i in range(n + 1):
        for j in range(i):
            c = aa[i] * bb[j] - bb[i] * aa[j]
            for k in range(i - j):
                B[i - 1 - k, j + k] += c
    return B


def synthesize_from_polynomials(p, h):
    """Entire synthesis up to rational Y uses coefficients only, no roots/eigs.

    Given p minimal polynomial of A=-D^2 on E(u), and q=D h(A)u,
    return F(lambda)=(N+t*p)/r with Y(s)=s F(-s^2) passive.
    """
    order = p.degree()
    r = P([0.0, 1.0]) * p.deriv() / order
    N = (r * h) % p
    H0 = -bezout(r, p)
    Bp = bezout(r, N)
    L = np.linalg.cholesky(H0)
    # M = L^-1 Bp L^-T; do not form inverse explicitly.
    tmp = np.linalg.solve(L, Bp)
    M = np.linalg.solve(L, tmp.T).T
    M = (M + M.T) / 2
    gersh = np.diag(M) + np.sum(np.abs(M), axis=1) - np.abs(np.diag(M))
    t = max(0.0, float(np.max(gersh)))
    return r, N, H0, Bp, M, t


def ftrim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def fmul(a, b):
    c = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return ftrim(c)


def frem(a, b):
    a = ftrim(a)
    while len(a) >= len(b):
        fac = a[-1] / b[-1]
        offset = len(a) - len(b)
        for j, x in enumerate(b):
            a[offset + j] -= fac * x
        a = ftrim(a)
    return a


def feval(a, x):
    out = F(0)
    for ai in reversed(a):
        out = out * x + ai
    return out


def exact_three_mode_check():
    # p=(lambda-1)(lambda-4)(lambda-9), target [1,-1,1].
    p = list(map(F, [-36, 49, -14, 1]))
    r = [F(0), F(49, 3), F(-28, 3), F(1)]
    h = [F(11, 5), F(-4, 3), F(2, 15)]
    N = frem(fmul(r, h), p)
    poles = [F(0), F(7, 3), F(7)]
    thresholds = [-feval(N, x) / feval(p, x) for x in poles]
    t = max([F(0)] + thresholds)
    rp = [r[j] * j for j in range(1, len(r))]
    weights = [-(feval(N, x) + t * feval(p, x)) / feval(rp, x) for x in poles]
    values = [(feval(N, x) + t * feval(p, x)) / feval(r, x) for x in map(F, [1, 4, 9])]
    print("EXACT H=3 (ascending coefficients)")
    print("p =", p)
    print("r =", r)
    print("h =", h)
    print("N =", N)
    print("pole thresholds =", thresholds)
    print("minimal t =", t)
    print("Foster positive coefficients [C, Gamma0, Gamma1, Gamma2] =", [t] + weights)
    print("matched values =", values)
    assert values == [F(1), F(-1), F(1)]
    assert min([t] + weights) >= 0


def randomized_checks():
    rng = np.random.default_rng(123)
    print("\nRANDOM CHECKS")
    for order in range(1, 7):
        nodes = np.linspace(0.2, 2, order) + rng.uniform(0, 0.05, order)
        targets = rng.normal(size=order)
        p = P.fromroots(nodes)
        h = P(np.linalg.solve(np.vander(nodes, order, increasing=True), targets))
        r, N, H0, Bp, M, t = synthesize_from_polynomials(p, h)
        # Independent verification below MAY use roots/eigenvalues.
        poles = r.roots()
        weights = -(N(poles) + t * p(poles)) / r.deriv()(poles)
        err = np.max(np.abs((N(nodes) + t * p(nodes)) / r(nodes) - targets))
        exact_t = max(0.0, float(eigh(Bp, H0, eigvals_only=True)[-1]))
        max_lmi = float(np.linalg.eigvalsh(Bp - t * H0)[-1])
        print(f"H={order}: t_bound={t:.9g}, t_min={exact_t:.9g}, "
              f"min_residue_weight={weights.min():.3g}, "
              f"interpolation_error={err:.3g}, max_LMI_eig={max_lmi:.3g}")
        assert t >= exact_t - 1e-7
        assert weights.min() >= -1e-7
        assert err <= 1e-7
        assert max_lmi <= 1e-7


def matrix_poly(poly, A):
    eye = np.eye(A.shape[0])
    out = np.zeros_like(A)
    for coeff in poly.coef[::-1]:
        out = out @ A + coeff * eye
    return out


def temporal_sample_check():
    """Identify p,h from time samples, then check synthesized current.

    Periodic trigonometric collocation D is written explicitly; no FFT.
    The known test current is constructed to have quadrature, in-phase,
    and generated components. Only (D,u,i) enter the reconstruction.
    """
    count = 25
    tt = 2 * np.pi * np.arange(count) / count
    jj = np.arange(count)
    delta = jj[:, None] - jj[None, :]
    nonzero = delta != 0
    D = np.zeros((count, count))
    D[nonzero] = 0.5 * (-1.0) ** delta[nonzero] / np.sin(np.pi * delta[nonzero] / count)
    A = -(D @ D)
    u = np.sin(tt) + 0.4 * np.cos(2 * tt) + 0.3 * np.sin(3 * tt)
    h_reference = P([11 / 5, -4 / 3, 2 / 15])
    q_reference = D @ matrix_poly(h_reference, A) @ u
    i = 0.8 * u + 0.1 * (A @ u) - q_reference + 0.2 * np.cos(4 * tt)
    # A rank test would identify H. This exact H=3 test uses known rank.
    H = 3
    E = np.column_stack([np.linalg.matrix_power(A, k) @ u for k in range(H)])
    coeff = np.linalg.lstsq(E, -(np.linalg.matrix_power(A, H) @ u), rcond=None)[0]
    p = P(np.r_[coeff, 1.0])
    K = D @ E
    h = P(np.linalg.lstsq(K, -i, rcond=None)[0])
    r, N, H0, Bp, M, t = synthesize_from_polynomials(p, h)
    # Evaluate rational response only on the excited cyclic subspace.
    Q, _ = np.linalg.qr(E)
    J = Q.T @ A @ Q
    q_synth = D @ Q @ np.linalg.solve(matrix_poly(r, J), matrix_poly(N + t * p, J) @ (Q.T @ u))
    q_projection = K @ h.coef
    err = np.linalg.norm(q_synth - q_projection) / np.linalg.norm(q_projection)
    referr = np.linalg.norm(q_projection - q_reference) / np.linalg.norm(q_reference)
    lossless = float(u @ q_synth / count)
    normal_error = np.linalg.norm(K.T @ (i + q_synth)) / (np.linalg.norm(K) * np.linalg.norm(i))
    print("\nTIME-SAMPLE CONSTRUCTION (no modal decomposition)")
    print("identified p coefficients =", p.coef)
    print("identified h coefficients =", h.coef)
    print("Gershgorin t =", t)
    print("rational current/projection relative error =", err)
    print("projection/reference relative error =", referr)
    print("mean u*q_synth =", lossless)
    print("normalized projection optimality residual =", normal_error)
    assert err < 1e-9
    assert referr < 1e-9
    assert abs(lossless) < 1e-9
    assert normal_error < 1e-9


if __name__ == "__main__":
    exact_three_mode_check()
    randomized_checks()
    temporal_sample_check()
