#!/usr/bin/env python3
"""Auditoria espectral exacta independiente de la busqueda temporal MATLAB.

Fuente ideal, soporte finito y LC serie positiva. Aisla todas las raices
estacionarias mediante Sturm y compara intervalos racionales de la mejora
de I_RMS^2 con C sola, L sola y los limites resonantes. No analiza la fuente
R-L acoplada. No importa ni ejecuta los certificados historicos.
"""
import ast
import csv
import json
import math
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# Reutiliza solo las funciones puras del certificado, sin ejecutar su campaña.
tree = ast.parse((ROOT/'analysis_incremental/quartic_lc_certificate.py').read_text())
names = {'trim', 'add', 'scale', 'mul', 'deriv', 'evaluate', 'primitive',
         'rational_remainder', 'sturm', 'variations'}
ns = dict(Fraction=F, reduce=__import__('functools').reduce,
          lcm=math.lcm, gcd=math.gcd)
exec(compile(ast.Module(body=[a for a in tree.body if isinstance(a, ast.FunctionDef)
                             and a.name in names], type_ignores=[]), '<exact-polynomials>', 'exec'), ns)
globals().update({name: ns[name] for name in names})


def interval_mul(a, b):
    values = [x*y for x in a for y in b]
    return min(values), max(values)


def interval_eval(coeff, x):
    out = (F(0), F(0))
    for v in reversed(coeff):
        lo, hi = interval_mul(out, x)
        out = (lo+v, hi+v)
    return out


def isolate_positive(poly):
    chain = sturm(poly)
    bound = F(1)+max(abs(F(v, poly[-1])) for v in poly[:-1])
    pending = [(F(0), bound)]
    roots = []
    while pending:
        lo, hi = pending.pop()
        count = variations(chain, lo)-variations(chain, hi)
        if not count:
            continue
        if count == 1 and hi-lo < F(1, 10**22):
            roots.append((lo, hi))
            continue
        mid = (lo+hi)/2
        assert evaluate(poly, mid) != 0, 'Raiz racional: aislar expresamente'
        pending += [(lo, mid), (mid, hi)]
    return sorted(roots)


def loads(case):
    if case == 'Jeltsema2014':
        n, u = [1, 5], [F(10), F(5)]
        g = [F(1, 1+4*k*k) for k in n]
        b = [F(-2*k, 1+4*k*k) for k in n]
        return n, u, g, b, F(1)
    if case == 'Czarnecki1987_ideal':
        n, u = [1, 5, 7], [F(100), F(7, 2), F(2)]
        return n, u, [F(2, 1+3*k*k) for k in n], [F(-2*k, 1+3*k*k) for k in n], F(3)
    n, u = [1, 2, 3], [F(100), F(50), F(30)]
    if case == 'Czarnecki1981_tabla':
        return n, u, list(map(F, ['.02749', '.01231', '.04264'])), list(map(F, ['-.07419', '.09060', '.26810'])), F(1)
    g, b = [], []
    for k in n:
        gp, bp = F(1, 1+36*k*k), F(k, 12)-F(6*k, 1+36*k*k)
        re, im = 1+(gp-k*bp)/2, (bp+k*gp)/2
        g.append((gp*re+bp*im)/(re*re+im*im))
        b.append((bp*re-gp*im)/(re*re+im*im))
    return n, u, g, b, F(1)


def audit(case):
    n, u, g, susceptance, factor = loads(case)
    p, a, b = [F(1)], [F(0)], [F(0)]
    for k in n:
        p = mul(p, [-k*k, 1])
    for k, uk, bk in zip(n, u, susceptance):
        quotient = [F(1)]
        for m in n:
            if m != k:
                quotient = mul(quotient, [-m*m, 1])
        a = add(a, scale(quotient, uk*uk*k*bk))
        b = add(b, scale(mul(quotient, quotient), uk*uk*k*k))
    stationary = primitive(add(scale(mul(deriv(a), b), 2), scale(mul(a, deriv(b)), -1)))
    base = sum(uk*uk*(gk*gk+factor*bk*bk) for uk, gk, bk in zip(u, g, susceptance))
    candidates = []
    for k, uk, bk in zip(n, u, susceptance):
        gain = factor*uk*uk*bk*bk
        candidates.append(dict(type='limite resonante', z=float(k*k), gain=(gain, gain)))
    for name, basis in [('C sola', n), ('L sola', [-F(1, k) for k in n])]:
        dot = sum(uk*uk*bk*h for uk, bk, h in zip(u, susceptance, basis))
        gain = factor*min(F(0), dot)**2/sum(uk*uk*h*h for uk, h in zip(u, basis))
        candidates.append(dict(type=name, gain=(gain, gain)))
    roots = isolate_positive(stationary)
    for bounds in roots:
        z = sum(bounds)/2
        aa, pp, bb = interval_eval(a, bounds), interval_eval(p, bounds), interval_eval(b, bounds)
        ap = interval_mul(aa, pp)
        assert ap[0]*ap[1] > 0 and bb[0] > 0
        if ap[0] > 0:  # alpha=-a*p/b debe ser positivo
            continue
        a2 = interval_mul(aa, aa)
        gain = (factor*a2[0]/bb[1], factor*a2[1]/bb[0])
        alpha = -float(evaluate(a, z)*evaluate(p, z)/evaluate(b, z))*math.sqrt(factor)
        candidates.append(dict(type='minimo interior candidato', z=float(z), gain=gain,
                               L_H=1/alpha, C_F=alpha/float(z),
                               root_interval=[str(x) for x in bounds]))
    candidates.sort(key=lambda x: sum(x['gain']), reverse=True)
    winner = candidates[0]
    assert all(winner['gain'][0] > x['gain'][1] for x in candidates[1:]), 'Maximo no separado'
    for row in candidates:
        lo, hi = row.pop('gain')
        row['gain_interval_exact'] = [str(lo), str(hi)]
        row['I_RMS_A'] = math.sqrt(float(base-(lo+hi)/2))
    return dict(case=case, passed=True, positive_stationary_roots=len(roots),
                stationary_coefficients_ascending=stationary,
                classification=('infimo no alcanzado; frontera resonante' if winner['type']=='limite resonante'
                                else 'minimo global LC positivo finito'),
                winner=winner, candidates=candidates,
                scope='Fuente ideal, soporte finito, LC serie positiva; auditoria espectral exacta independiente.')


def main():
    records = [audit(c) for c in ['Jeltsema2014', 'Czarnecki1987_ideal', 'Czarnecki1981', 'Czarnecki1981_tabla']]
    folder = ROOT/'paper/comparacion'
    result = dict(passed=True, records=records,
                  proof='Sturm con racionales: todas las raices positivas. Comparacion por intervalos racionales del beneficio con todas las fronteras.')
    (folder/'fronteras_exactas.json').write_text(json.dumps(result, indent=2)+'\n')
    with (folder/'fronteras_exactas.csv').open('w') as stream:
        writer = csv.writer(stream)
        writer.writerow(['Caso', 'Clasificacion', 'Is_limite_o_minimo_A', 'z', 'L_H', 'C_F'])
        for r in records:
            w = r['winner']
            writer.writerow([r['case'], r['classification'], w['I_RMS_A'], w.get('z'), w.get('L_H'), w.get('C_F')])
            print(r['case'], r['classification'], w['I_RMS_A'], w.get('z'))


if __name__ == '__main__':
    main()
