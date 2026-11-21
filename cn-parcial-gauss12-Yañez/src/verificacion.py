"""Verificación independiente"""

import math
from fractions import Fraction

import mpmath

from gauss import resolver_con_LU, resolver_transpuesta_con_LU

mpmath.mp.dps = 50


def resolver_exacto(A_fr, b_fr):
    """Eliminación gaussiana en aritmética racional exacta (referencia, solo para verificar)."""
    n = len(A_fr)
    M = [fila[:] + [b_fr[i]] for i, fila in enumerate(A_fr)]
    for k in range(n):
        p = next(i for i in range(k, n) if M[i][k] != 0)
        M[k], M[p] = M[p], M[k]
        for i in range(k + 1, n):
            m = M[i][k] / M[k][k]
            M[i] = [M[i][j] - m * M[k][j] for j in range(n + 1)]
    x = [Fraction(0)] * n
    for i in range(n - 1, -1, -1):
        x[i] = (M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))) / M[i][i]
    return x


def residuo_exacto(A_fr, b_fr, x):
    """r = b − A·x con x tomado tal cual (float o Fraction) y A, b racionales exactos: sin redondeo."""
    xf = [Fraction(v) for v in x]
    return [b_fr[i] - sum(A_fr[i][j] * xf[j] for j in range(len(xf))) for i in range(len(xf))]


def norma_inf(v):
    return max(abs(t) for t in v)


def a_mp(fr):
    return mpmath.mpf(fr.numerator) / mpmath.mpf(fr.denominator)


def fmt_mp(fr, cifras=6):
    return mpmath.nstr(a_mp(Fraction(fr)), cifras)


def digitos_correctos(x, xref):
    """Dígitos significativos correctos por componente: floor(−log10 |x−xref|/|xref|)."""
    out = []
    for xi, ri in zip(x, xref):
        err = abs(Fraction(xi) - ri)
        if err == 0:
            out.append(17)
            continue
        rel = float(err / abs(ri))
        out.append(max(0, min(17, math.floor(-math.log10(rel)))))
    return out


def error_relativo_norma(x, xref):
    num = max(abs(Fraction(a) - b) for a, b in zip(x, xref))
    den = max(abs(b) for b in xref)
    return float(num / den)


def refinar(f, A_fr, b_fr, x0, iteraciones=3, extendido=False):
    """Refinamiento iterativo reutilizando L y U (sin refactorizar).

    Residuo exacto (racional) y corrección d = (LU)⁻¹ r en doble. Si extendido=False el iterado se guarda
    en doble; si True se acumula exactamente (Fraction) para ver la convergencia bajo ε.
    """
    x = [Fraction(v) for v in x0] if extendido else list(x0)
    hist = []
    for it in range(1, iteraciones + 1):
        r = residuo_exacto(A_fr, b_fr, x)
        d = resolver_con_LU(f, [float(t) for t in r])
        x = [Fraction(a) + Fraction(b) for a, b in zip(x, d)] if extendido else [a + b for a, b in zip(x, d)]
        hist.append({"it": it, "res_inf": norma_inf(r), "d": d, "x": list(x)})
    return x, hist


#condicionamiento

def norma1(A):
    n = len(A)
    return max(sum(abs(A[i][j]) for i in range(n)) for j in range(n))


def estimar_norma1_inversa(f, n, max_iter=5):
    """Estimador de Hager (1984) de ‖A⁻¹‖₁ usando solo resoluciones con L y U (nunca se forma A⁻¹)."""
    x = [1.0 / n] * n
    est = 0.0
    for _ in range(max_iter):
        y = resolver_con_LU(f, x)
        est = sum(abs(v) for v in y)
        xi = [1.0 if v >= 0 else -1.0 for v in y]
        z = resolver_transpuesta_con_LU(f, xi)
        j = max(range(n), key=lambda t: abs(z[t]))
        if abs(z[j]) <= sum(zi * xi_ for zi, xi_ in zip(z, x)):
            break
        x = [0.0] * n
        x[j] = 1.0
    return est


def kappa1_estimado(A, f):
    return norma1(A) * estimar_norma1_inversa(f, len(A))


def casi_dependencia(A, i, j, k):
    """Ajusta la fila i como α·fila_j + β·fila_k (mínimos cuadrados 2x2, ecuaciones normales).

    Devuelve (α, β, residuo relativo en norma 2). Un residuo pequeño delata casi dependencia lineal.
    """
    def dot(u, v):
        return sum(a * b for a, b in zip(u, v))

    u, v, w = A[j], A[k], A[i]
    a11, a12, a22 = dot(u, u), dot(u, v), dot(v, v)
    y1, y2 = dot(u, w), dot(v, w)
    d = a11 * a22 - a12 * a12
    alfa = (y1 * a22 - y2 * a12) / d
    beta = (a11 * y2 - a12 * y1) / d
    res = [w[t] - alfa * u[t] - beta * v[t] for t in range(len(w))]
    return alfa, beta, math.sqrt(dot(res, res)) / math.sqrt(dot(w, w))
