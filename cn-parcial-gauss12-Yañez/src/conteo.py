"""Instrumentación y reporte del costo operacional (los números salen de gauss.Contador)."""

REF_MINIMO = 1486          # cota de referencia de la prueba: eliminación densa bien implementada
REF_PIVOTEO_EXTRA = 66     # comparaciones de pivoteo parcial según la prueba
REF_RECALCULADO = 1750     # multiplicadores recalculados por columna (según la prueba)
REF_INVERSA = 7000         # vía matriz inversa (según la prueba)


def formula_de_control(n):
    """Fórmula cerrada SOLO como contraste del contador (no es de donde sale el conteo)."""
    mul_div = sum((n - k) * (n - k + 2) for k in range(1, n)) + (n - 1) * n // 2 + n
    sumas = sum((n - k) * (n - k + 1) for k in range(1, n)) + (n - 1) * n // 2
    return mul_div, sumas


def tabla_pasos_md(c):
    n = c.n
    filas = ["| Paso k | Divisiones (multiplicadores) | Multiplicaciones | Sumas/restas | Total flops |",
             "| ---: | ---: | ---: | ---: | ---: |"]
    for k, p in enumerate(c.elim, start=1):
        filas.append(f"| {k} | {p['div']} | {p['mul']} | {p['add']} | {p['div'] + p['mul'] + p['add']} |")
    e = c.total_elim()
    sd = sum(p["div"] for p in c.elim); sm = sum(p["mul"] for p in c.elim); sa = sum(p["add"] for p in c.elim)
    filas.append(f"| **Eliminación (11 pasos)** | **{sd}** | **{sm}** | **{sa}** | **{e}** |")
    r = c.retro
    filas.append(f"| Sustitución regresiva | {r['div']} | {r['mul']} | {r['add']} | {c.total_retro()} |")
    filas.append(f"| **TOTAL (sin pivoteo)** | **{sd + r['div']}** | **{sm + r['mul']}** | **{sa + r['add']}** | **{c.total_base()}** |")
    return "\n".join(filas)


def tabla_pivoteo_md(c):
    filas = ["| Paso k | Divisiones de razón (escalado) | Comparaciones | Intercambios de fila |",
             "| ---: | ---: | ---: | ---: |"]
    for k, p in enumerate(c.piv_paso, start=1):
        filas.append(f"| {k} | {p['div']} | {p['cmp']} | {p['swaps']} |")
    filas.append(f"| Escalas s_i (máx por fila) | – | {c.piv_escala_cmp} | – |")
    filas.append(f"| **Sobrecosto total de estabilidad** | **{sum(p['div'] for p in c.piv_paso)}** | "
                 f"**{c.piv_escala_cmp + sum(p['cmp'] for p in c.piv_paso)}** | **{c.total_swaps()}** |")
    return "\n".join(filas)


def _desv(valor, ref):
    return f"{100.0 * (valor - ref) / ref:+.1f} %"


def comparar_md(base, con_pivoteo, recalculado):
    """base: total sin pivoteo; con_pivoteo: total + sobrecosto (escalado); recalculado: variante redundante."""
    ref_p = REF_MINIMO + REF_PIVOTEO_EXTRA
    return "\n".join([
        "| Estrategia | Flops medidos | Cota de la prueba | Desviación | ¿Dentro de ±10 %? |",
        "| --- | ---: | ---: | ---: | :---: |",
        f"| Eliminación densa, multiplicadores reutilizados | {base} | ≈ {REF_MINIMO} | {_desv(base, REF_MINIMO)} | {'sí' if abs(base-REF_MINIMO)<=0.1*REF_MINIMO else 'no'} |",
        f"| Lo mismo + pivoteo (búsqueda y escalas) | {con_pivoteo} | ≈ {ref_p} (+66 cmp) | {_desv(con_pivoteo, ref_p)} | {'sí' if abs(con_pivoteo-ref_p)<=0.1*ref_p else 'no'} |",
        f"| Multiplicadores recalculados por columna (variante penalizada) | {recalculado} | ≈ {REF_RECALCULADO} | {_desv(recalculado, REF_RECALCULADO)} | exceso: +{recalculado-base} flops sobre el mínimo |",
        f"| Vía matriz inversa (no implementada, descalifica C2) | – | ≈ {REF_INVERSA} o más | – | – |",
    ])


def reporte_md(n, conts):
    """conts: dict con 'ninguno', 'parcial', 'escalado', 'recalculado' -> Contador."""
    base = conts["ninguno"].total_base()
    assert all(conts[m].total_base() == base for m in ("parcial", "escalado")), "el pivoteo no debe cambiar el costo base"
    md, sm = formula_de_control(n)
    assert md + sm == base, "el contador no coincide con la fórmula de control"
    esc = conts["escalado"]; par = conts["parcial"]; rec = conts["recalculado"]
    partes = [
        "# Conteo de operaciones (n = 12)",
        "",
        "**Convención.** 1 flop = 1 operación aritmética elemental (+, −, ×, ÷) en punto flotante. "
        "Los números salen de `Contador`, que se incrementa dentro de los bucles de `src/gauss.py`. "
        "El elemento a(i,k) que se anula no se calcula (se fija en 0); cada multiplicador se calcula una sola vez "
        "y se reutiliza en toda la fila y en b. `abs()` y los intercambios de fila no se cuentan como flops "
        "(los intercambios se informan aparte). Si además se calculase explícitamente a(i,k) − m·a(k,k) "
        "(66 multiplicaciones y 66 restas más) el total sería 1486, que es la cota de referencia de la prueba.",
        "",
        "## Eliminación + sustitución regresiva (costo base, idéntico con o sin pivoteo)",
        "",
        tabla_pasos_md(esc),
        "",
        f"Control cruzado: la fórmula cerrada da {md} multiplicaciones/divisiones y {sm} sumas/restas "
        f"({md + sm} en total) y coincide con el contador.",
        "",
        "## Sobrecosto de estabilidad: pivoteo parcial escalado (contado aparte)",
        "",
        tabla_pivoteo_md(esc),
        "",
        f"Pivoteo parcial simple (sin escalar): {sum(p['cmp'] for p in par.piv_paso)} comparaciones, "
        f"{par.total_swaps()} intercambios, sin divisiones ni escalas.",
        "",
        "## Comparación con las cotas de la prueba",
        "",
        comparar_md(base, base + esc.sobrecosto_pivoteo(), rec.total_base()),
        "",
    ]
    return "\n".join(partes)
