"""Punto de entrada único: reproduce todo el contenido de resultados/ con `python src/main.py`."""
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import conteo
import gauss
import verificacion as ver

RAIZ = Path(__file__).resolve().parent.parent

EPS = 2.0 ** -52          # ε de máquina en doble precisión (≈ 2.22e-16)

MODOS = ("ninguno", "parcial", "escalado")


def escribir_traza(ruta, corridas):
    with open(ruta, "w", encoding="utf-8") as f:
        for modo in ("escalado", "ninguno"):
            fz = corridas[modo]
            f.write("=" * 100 + "\n")
            f.write(f"TRAZA DE ELIMINACIÓN — modo: {'pivoteo parcial escalado' if modo == 'escalado' else 'sin pivoteo'}\n")
            f.write("Matriz aumentada [A | b] tras cada paso (9 cifras significativas en la impresión; "
                    "los cálculos usan doble precisión completa).\n")
            f.write("=" * 100 + "\n\n")
            for k, (snap, piv) in enumerate(zip(fz.traza, fz.pivotes), start=1):
                f.write(f"--- Paso {k}: pivote a({k},{k}) = {piv[2]:+.9e} (fila original {piv[1]})\n")
                for fila in snap:
                    f.write(" ".join(f"{v:+.9e}" for v in fila) + "\n")
                bajo = sum(1 for i in range(len(snap)) for j in range(min(i, k)) if snap[i][j] != 0.0)
                f.write(f"    elementos no nulos bajo la diagonal en las columnas 1..{k}: {bajo}\n\n")
            n = len(fz.U)
            no_nulos = sum(1 for i in range(n) for j in range(i) if fz.U[i][j] != 0.0)
            f.write(f"RESULTADO: elementos no nulos bajo la diagonal = {no_nulos} de 66 (debe ser 0)\n\n")
            assert no_nulos == 0


def main():
    A_txt = gauss.leer_matriz_texto(RAIZ / "datos" / "matriz_a.csv")
    b_txt = gauss.leer_vector_texto(RAIZ / "datos" / "vector_b.csv")
    n = len(A_txt)
    assert n == 12 and all(len(r) == 12 for r in A_txt) and len(b_txt) == 12
    A, b = gauss.a_float(A_txt), [float(s) for s in b_txt]
    A_fr, b_fr = gauss.a_fraccion(A_txt), [Fraction(s) for s in b_txt]

    # referencia exacta (racional) y las corridas del mismo código

    xref = ver.resolver_exacto(A_fr, b_fr)
    corridas = {m: gauss.eliminar(A, b, m) for m in MODOS}
    recalc = gauss.eliminar(A, b, "ninguno", reusar_multiplicador=False)
    conts = {m: corridas[m].contador for m in MODOS}
    conts["recalculado"] = recalc.contador

    #traza y conteo

    escribir_traza(RAIZ / "resultados" / "traza_eliminacion.txt", corridas)
    (RAIZ / "resultados" / "conteo_flops.md").write_text(conteo.reporte_md(n, conts), encoding="utf-8")

    #comparación de corridas

    f = corridas["escalado"]
    x0 = f.x
    lineas = []
    P = lineas.append
    P("RESIDUO Y VERIFICACIÓN — todos los residuos usan A y b exactos (racionales) y x tal cual (sin redondeo)\n")
    P("== T2. Comparación de las corridas del mismo código ==")
    P(f"{'modo':<10} {'mult. máx':>12} {'bits mult':>10} {'crecimiento':>12} {'‖r‖∞':>12} {'err. rel. (norma ∞)':>20} {'díg. min':>9} {'bits perdidos':>14}")
    resumen = {}
    for m in MODOS:
        fz = corridas[m]
        r = ver.residuo_exacto(A_fr, b_fr, fz.x)
        err = ver.error_relativo_norma(fz.x, xref)
        dig = ver.digitos_correctos(fz.x, xref)
        perdidos = math.log2(err / 2.0 ** -53) if err > 0 else 0.0
        resumen[m] = (err, dig)
        P(f"{m:<10} {fz.mult_max:>12.4e} {math.log2(fz.mult_max):>10.1f} {fz.crecimiento:>12.4e} "
          f"{ver.fmt_mp(ver.norma_inf(r), 4):>12} {err:>20.4e} {min(dig):>9d} {perdidos:>14.1f}")
    P("bits mult = log2(multiplicador máximo): bits de mantisa del operando pequeño que el alineamiento de exponentes puede empujar fuera de los 53 bits.")
    P("bits perdidos = log2(error relativo / 2^-53): pérdida efectiva medida frente a la referencia exacta.\n")

    #residuo del método principal

    r0 = ver.residuo_exacto(A_fr, b_fr, x0)
    P("== T4. Residuo r = b − A·x (pivoteo parcial escalado), precisión extendida (racional exacto) ==")
    for i, ri in enumerate(r0, start=1):
        P(f"r[{i:2d}] = {ver.fmt_mp(ri, 10):>18}")
    P(f"‖r‖∞ = {ver.fmt_mp(ver.norma_inf(r0), 10)}")
    P(f"‖b‖∞ = {ver.fmt_mp(max(abs(t) for t in b_fr), 10)}")
    normA = max(sum(abs(t) for t in fila) for fila in A_fr)
    normx = max(abs(Fraction(v)) for v in x0)
    P(f"‖A‖∞ = {ver.fmt_mp(normA, 6)}   ‖x‖∞ = {ver.fmt_mp(normx, 6)}")
    P(f"error inverso relativo ‖r‖∞ / (‖A‖∞·‖x‖∞) = {ver.fmt_mp(ver.norma_inf(r0) / (normA * normx), 4)}   (u = 2^-53 ≈ 1.11e-16)\n")

    #refinamiento iterativo

    P("== Refinamiento iterativo (reutiliza L y U; residuo exacto, corrección en doble) ==")
    xr, hist = ver.refinar(f, A_fr, b_fr, x0, iteraciones=3, extendido=False)
    xe, hist_e = ver.refinar(f, A_fr, b_fr, x0, iteraciones=3, extendido=True)
    P(f"{'iter':<6} {'‖r‖∞ (antes)':>16} {'‖d‖∞/‖x‖∞':>14} {'err. rel. iterado en doble':>28} {'err. rel. iterado extendido':>29}")
    P(f"{0:<6} {ver.fmt_mp(ver.norma_inf(r0), 4):>16} {'–':>14} {ver.error_relativo_norma(x0, xref):>28.4e} {ver.error_relativo_norma(x0, xref):>29.4e}")
    xnorm = max(abs(v) for v in x0)
    for h, he in zip(hist, hist_e):
        dn = max(abs(v) for v in h["d"]) / xnorm
        P(f"{h['it']:<6} {ver.fmt_mp(h['res_inf'], 4):>16} {dn:>14.4e} {ver.error_relativo_norma(h['x'], xref):>28.4e} {ver.error_relativo_norma(he['x'], xref):>29.4e}")
    P("")

    #condicionamiento

    k1 = ver.kappa1_estimado(A, f)
    k2 = float(np.linalg.cond(np.array(A), 2))      # control (SVD), no resuelve el sistema
    cota = k2 * EPS
    P("== Condicionamiento ==")
    P(f"κ₁(A) estimado (Hager sobre L y U, sin formar A⁻¹) = {k1:.4e}")
    P(f"κ₂(A) control con numpy.linalg.cond            = {k2:.4e}")
    P(f"ε (doble precisión, 2^-52)                     = {EPS:.4e}")
    P(f"cota κ₂·ε sobre ‖δx‖/‖x‖                       = {cota:.4e}  -> {max(0, -math.log10(cota)):.1f} dígitos garantizados por la cota")
    P(f"error real (escalado, norma ∞)                 = {resumen['escalado'][0]:.4e}  -> {-math.log10(resumen['escalado'][0]):.1f} dígitos reales")
    escalas = [max(abs(v) for v in fila) for fila in A]
    Aeq = np.array([[v / escalas[i] for v in fila] for i, fila in enumerate(A)])
    k2eq = float(np.linalg.cond(Aeq, 2))            # control (SVD) tras equilibrar filas
    P(f"κ₂ tras equilibrar filas (cada fila / su máx |a_ij|) = {k2eq:.4e}  -> κ_eq·ε = {k2eq * EPS:.4e} ({-math.log10(k2eq * EPS):.1f} dígitos)")
    alfa, beta, rel = ver.casi_dependencia(A, 6, 2, 4)
    P(f"casi dependencia: fila 7 ≈ {alfa:.6g}·fila 3 + {beta:.6g}·fila 5, residuo relativo (norma 2) = {rel:.3e}")
    P("")

    #control con biblioteca

    xn = np.linalg.solve(np.array(A), np.array(b)).tolist()
    P("== Control con numpy.linalg.solve (no sustituye la implementación propia) ==")
    P(f"error relativo numpy vs exacta (norma ∞)       = {ver.error_relativo_norma(xn, xref):.4e}  (dígitos mínimos correctos: {min(ver.digitos_correctos(xn, xref))})")
    P(f"máx |x_propia − x_numpy| / ‖x‖∞                = {max(abs(a - c) for a, c in zip(x0, xn)) / xnorm:.4e}")
    (RAIZ / "resultados" / "residuo.txt").write_text("\n".join(lineas) + "\n", encoding="utf-8")

    #solución con dígitos defendibles
    
    d1 = hist[0]["d"]
    dig_ver = ver.digitos_correctos(x0, xref)
    filas = ["i,x,digitos_afirmados,digitos_verificados_vs_exacta,x_refinada"]
    for i in range(n):
        est = abs(d1[i]) / abs(x0[i]) if d1[i] != 0 else 0.0
        afirmados = 15 if est == 0 else max(1, min(15, math.floor(-math.log10(est))))
        filas.append(f"{i + 1},{x0[i]:.{afirmados - 1}e},{afirmados},{dig_ver[i]},{xr[i]:.16e}")
    (RAIZ / "resultados" / "solucion.csv").write_text("\n".join(filas) + "\n", encoding="utf-8")

    print("\n".join(lineas))
    print("\nsolucion.csv:")
    print("\n".join(filas))
    print(f"\nFlops base: {conts['escalado'].total_base()}  sobrecosto pivoteo escalado: {conts['escalado'].sobrecosto_pivoteo()}")


if __name__ == "__main__":
    main()
