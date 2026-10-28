"""Cauer-I LC synthesis using polynomial divisions only.

Inputs are numerator and denominator coefficients in ASCENDING powers.
There is no root calculation, eigendecomposition, or Fourier transform.
The routines are intended for normalized, low-order polynomial data whose
lossless positive-real property has already been established. Float failures
raise an error rather than silently returning negative/non-LC components.

Exact Fraction inputs are supported for algebraic verification. The output
alternates shunt C and series L, beginning at the input port. A final shunt C
ends at the return conductor; a final series L is terminated by the return
conductor. This implementation expects strictly positive initial C; increase
the free null-network parameter tau when needed to ensure that condition.
"""
from fractions import Fraction
import math


def _is_exact(values):
    return all(isinstance(x, (int, Fraction)) for x in values)


def _trim(a, exact, rtol, atol):
    a = list(a)
    if not a:
        return [Fraction(0) if exact else 0.0]
    scale = max(abs(float(x)) for x in a)
    threshold = 0.0 if exact else atol + rtol * scale
    while len(a) > 1 and (a[-1] == 0 if exact else abs(a[-1]) <= threshold):
        a.pop()
    return a


def _zero(a, exact, scale, rtol, atol):
    if exact:
        return all(x == 0 for x in a)
    return max(abs(x) for x in a) <= atol + rtol * scale


def rational_y_from_effective(n_effective, d_effective):
    """Build Y(s)=s*n_effective(-s²)/d_effective(-s²).

    Common factors at s=0 are cancelled directly, without factorization.
    Returns ascending coefficients (numerator, denominator).
    """
    n_effective = list(n_effective)
    d_effective = list(d_effective)
    exact = _is_exact(n_effective + d_effective)
    zero = Fraction(0) if exact else 0.0
    n = [zero] * (2 * len(n_effective))
    d = [zero] * (2 * len(d_effective) - 1)
    for k, coeff in enumerate(n_effective):
        n[2 * k + 1] = (-1) ** k * coeff
    for k, coeff in enumerate(d_effective):
        d[2 * k] = (-1) ** k * coeff
    while len(n) > 1 and len(d) > 1 and n[0] == 0 and d[0] == 0:
        n.pop(0)
        d.pop(0)
    return n, d


def cauer_lc(num, den, rtol=1e-11, atol=1e-13, max_elements=None):
    """Extract a positive Cauer-I ladder from the rational admittance num/den.

    Output is a list of {'connection': 'shunt'|'series',
                         'element': 'C'|'L', 'value': scalar}.
    The value is capacitance for shunt C or inductance for series L.
    With Fraction coefficients, all calculations and returned values are exact.

    Assumes the first nonzero quotient is proportional to s. The conservative
    tau from the Bezout construction can be increased slightly to ensure a
    strictly positive first capacitor and avoid exact pole cancellations.
    """
    exact = _is_exact(list(num) + list(den))
    if exact:
        n = [Fraction(x) for x in num]
        d = [Fraction(x) for x in den]
    else:
        n = [float(x) for x in num]
        d = [float(x) for x in den]
    n = _trim(n, exact, rtol, atol)
    d = _trim(d, exact, rtol, atol)
    if max_elements is None:
        max_elements = len(n) + len(d) + 2
    ladder = []
    for index in range(max_elements):
        if _zero(d, exact, 1.0, rtol, atol):
            raise ValueError("Zero denominator during Cauer synthesis.")
        if len(n) != len(d) + 1:
            raise ValueError(
                f"Expected numerator degree one above denominator at step {index}; "
                f"got {len(n)-1} and {len(d)-1}. Check cancellations/scaling."
            )
        coefficient = n[-1] / d[-1]
        if coefficient <= 0 or not math.isfinite(float(coefficient)):
            raise ValueError(f"Nonpositive or invalid LC coefficient at step {index}: {coefficient}")
        ladder.append({
            "connection": "shunt" if index % 2 == 0 else "series",
            "element": "C" if index % 2 == 0 else "L",
            "value": coefficient,
        })
        # Polynomial long division: n/d = coefficient*s + remainder/d.
        remainder = list(n)
        for j, dj in enumerate(d):
            remainder[j + 1] -= coefficient * dj
        scale = max(max(abs(float(x)) for x in n),
                    max(abs(float(coefficient * x)) for x in d))
        # Leading term cancels by construction. Do not let rounding create an
        # artificial high-degree residual; require it to be negligible first.
        if exact:
            if remainder[-1] != 0:
                raise ArithmeticError("Exact polynomial cancellation failed.")
        elif abs(remainder[-1]) > atol + rtol * scale:
            raise ArithmeticError("Leading polynomial cancellation lost precision.")
        remainder.pop()
        if _zero(remainder, exact, scale, rtol, atol):
            return ladder
        remainder = _trim(remainder, exact, rtol, atol)
        # Invert the residual: next stage is impedance, then admittance, etc.
        n, d = d, remainder
    raise ValueError("Cauer synthesis exceeded the expected finite order.")


def evaluate_ladder(ladder, s):
    """Evaluate the input admittance by the continued fraction, from its end."""
    if not ladder:
        return 0.0 * s
    result = ladder[-1]["value"] * s
    for element in reversed(ladder[:-1]):
        result = element["value"] * s + 1 / result
    return result


def evaluate_polynomial(a, s):
    result = 0.0 * s
    for coeff in reversed(a):
        result = result * s + coeff
    return result


def demo():
    # Exact normalized example from verify_bezout_closure.py, H=3, tau=3.
    # p=lambda³-14lambda²+49lambda-36, N=2lambda²-6lambda+12.
    # r=lambda³-(28/3)lambda²+(49/3)lambda.
    # N + 3p = 3lambda³-40lambda²+141lambda-96.
    effective_n = [Fraction(-96), Fraction(141), Fraction(-40), Fraction(3)]
    effective_d = [Fraction(0), Fraction(49, 3), Fraction(-28, 3), Fraction(1)]
    num, den = rational_y_from_effective(effective_n, effective_d)
    ladder = cauer_lc(num, den)
    expected = [Fraction(3), Fraction(1, 12), Fraction(36, 5),
                Fraction(5, 96), Fraction(48, 5), Fraction(5, 144)]
    assert [x["value"] for x in ladder] == expected
    print("Exact Cauer-I ladder (normalized units):")
    for j, item in enumerate(ladder, 1):
        print(j, item["connection"], item["element"], item["value"])
    points = [0.1, 0.7, 1 + 1j, 0.2 + 2j, 0.5j, 1j, 2j, 3j, 4j]
    errors = []
    for s in points:
        rational = evaluate_polynomial(num, s) / evaluate_polynomial(den, s)
        reconstructed = evaluate_ladder(ladder, s)
        errors.append(abs(reconstructed - rational) / max(1.0, abs(rational)))
    print("Max relative reconstruction error:", max(errors))
    assert max(errors) < 1e-12
    # Re-run in ordinary floating-point arithmetic for notebook integration.
    float_ladder = cauer_lc(list(map(float, num)), list(map(float, den)))
    float_errors = []
    for s in points:
        rational = evaluate_polynomial(num, s) / evaluate_polynomial(den, s)
        reconstructed = evaluate_ladder(float_ladder, s)
        float_errors.append(abs(reconstructed - rational) / max(1.0, abs(rational)))
    print("Max float reconstruction error:", max(float_errors))
    assert max(float_errors) < 1e-11
    for j in [1, 2, 3]:
        print(f"Y(j*{j}) =", evaluate_ladder(ladder, 1j * j))


if __name__ == "__main__":
    demo()
