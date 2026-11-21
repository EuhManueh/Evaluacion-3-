"""Eliminación gaussiana densa con y sin pivoteo, instrumentada con contador de flops.

Convención de conteo (se declara también en informe.md y resultados/conteo_flops.md):
  * 1 flop = 1 operación aritmética elemental (+, -, *, /) sobre números de punto flotante.
  * abs() y los intercambios de filas no son flops; las comparaciones de la búsqueda de
    pivote y las divisiones del pivoteo escalado se cuentan APARTE (sobrecosto de estabilidad).
  * El elemento a(i,k) que se anula en el paso k NO se calcula: se fija en 0 por construcción
    (ahorro de n(n-1)/2 multiplicaciones y restas). Cada multiplicador m(i,k) = a(i,k)/a(k,k)
    se calcula UNA sola vez y se reutiliza en toda la fila (columnas k+1..n y lado derecho b).
  * Las actualizaciones de b se hacen dentro de la eliminación (columna aumentada).

El código NO usa ninguna rutina de biblioteca para resolver el sistema.
"""
from fractions import Fraction


#datos
def leer_matriz_texto(ruta):
    """Devuelve las 12 filas como listas de cadenas tal cual están en el CSV (sin perder cifras)."""
    with open(ruta, encoding="utf-8") as f:
        return [linea.strip().split(",") for linea in f if linea.strip()]


def leer_vector_texto(ruta):
    with open(ruta, encoding="utf-8") as f:
        return [linea.strip() for linea in f if linea.strip()]


def a_float(matriz_txt):
    return [[float(s) for s in fila] for fila in matriz_txt]


def a_fraccion(matriz_txt):
    """Valores decimales exactos como racionales (para residuo y referencia exacta)."""
    return [[Fraction(s) for s in fila] for fila in matriz_txt]



# contador

class Contador:
    """Contador de operaciones. Se incrementa dentro de los bucles de gauss.py (no es una fórmula)."""

    def __init__(self, n):
        self.n = n
        self.elim = [{"div": 0, "mul": 0, "add": 0} for _ in range(n - 1)]     # por paso k
        self.retro = {"div": 0, "mul": 0, "add": 0}                            # sustitución regresiva
        self.piv_paso = [{"cmp": 0, "div": 0, "swaps": 0} for _ in range(n - 1)]
        self.piv_escala_cmp = 0   # comparaciones para hallar la escala s_i de cada fila

    # totales
    def total_elim(self):
        return sum(p["div"] + p["mul"] + p["add"] for p in self.elim)

    def total_retro(self):
        return self.retro["div"] + self.retro["mul"] + self.retro["add"]

    def total_base(self):
        """Costo de eliminación + sustitución (sin sobrecosto de pivoteo)."""
        return self.total_elim() + self.total_retro()

    def sobrecosto_pivoteo(self):
        return (self.piv_escala_cmp
                + sum(p["cmp"] + p["div"] for p in self.piv_paso))

    def total_swaps(self):
        return sum(p["swaps"] for p in self.piv_paso)


#eliminación

class Factorizacion:
    """Resultado de P·A = L·U más lo necesario para reutilizar L y U con otros lados derechos."""

    def __init__(self, U, L, perm, c, x, mult_max, crecimiento, traza, pivotes, contador, modo):
        self.U, self.L, self.perm, self.c, self.x = U, L, perm, c, x
        self.mult_max, self.crecimiento = mult_max, crecimiento
        self.traza, self.pivotes, self.contador, self.modo = traza, pivotes, contador, modo


def _max_abs(M):
    return max(abs(v) for fila in M for v in fila)


def eliminar(A, b, modo="ninguno", reusar_multiplicador=True):
    """Triangulariza [A|b] y resuelve por sustitución regresiva.

    modo: 'ninguno' | 'parcial' (máx |a_ik|) | 'escalado' (máx |a_ik|/s_i, s_i = máx |a_ij| original).
    reusar_multiplicador=False recalcula m(i,k) por columna (variante redundante, solo para medir
    el exceso de costo que penaliza la sección T1).
    """
    n = len(A)
    M = [fila[:] for fila in A]
    rhs = b[:]
    L = [[0.0] * n for _ in range(n)]
    perm = list(range(n))
    cont = Contador(n)

    escala = None
    if modo == "escalado":
        escala = []
        for i in range(n):
            s = abs(M[i][0])
            for j in range(1, n):
                cont.piv_escala_cmp += 1          # una comparación por elemento tras el primero
                if abs(M[i][j]) > s:
                    s = abs(M[i][j])
            escala.append(s)

    max0 = _max_abs(M)
    mult_max = 0.0
    crecimiento = 1.0
    traza, pivotes = [], []

    for k in range(n - 1):

        #búsqueda de pivote (sobrecosto de estabilidad, contado aparte)

        p = k
        if modo != "ninguno":
            if modo == "escalado":
                cont.piv_paso[k]["div"] += 1      # razón del candidato inicial
                mejor = abs(M[k][k]) / escala[k]
            else:
                mejor = abs(M[k][k])
            for i in range(k + 1, n):
                if modo == "escalado":
                    cont.piv_paso[k]["div"] += 1
                    r = abs(M[i][k]) / escala[i]
                else:
                    r = abs(M[i][k])
                cont.piv_paso[k]["cmp"] += 1
                if r > mejor:
                    mejor, p = r, i
            if p != k:
                M[k], M[p] = M[p], M[k]
                rhs[k], rhs[p] = rhs[p], rhs[k]
                perm[k], perm[p] = perm[p], perm[k]
                if escala is not None:
                    escala[k], escala[p] = escala[p], escala[k]

                for j in range(k):                # multiplicadores ya guardados viajan con la fila
                    L[k][j], L[p][j] = L[p][j], L[k][j]
                cont.piv_paso[k]["swaps"] += 1
        pivotes.append((k + 1, perm[k] + 1, M[k][k]))


        # eliminación del paso k: sólo columnas k+1..n-1 y b (la columna k queda en 0)

        pk = cont.elim[k]
        for i in range(k + 1, n):
            m = M[i][k] / M[k][k]
            pk["div"] += 1
            L[i][k] = m
            if abs(m) > mult_max:
                mult_max = abs(m)
            for j in range(k + 1, n):
                if not reusar_multiplicador:

                    m = M[i][k] / M[k][k]         # recálculo redundante (variante penalizada)
                    pk["div"] += 1
                M[i][j] = M[i][j] - m * M[k][j]
                pk["mul"] += 1
                pk["add"] += 1
            if not reusar_multiplicador:
                m = M[i][k] / M[k][k]
                pk["div"] += 1
            rhs[i] = rhs[i] - m * rhs[k]
            pk["mul"] += 1
            pk["add"] += 1

            M[i][k] = 0.0                         # anulado por construcción, sin flops

        crecimiento = max(crecimiento, _max_abs(M) / max0)   # análisis, no cuenta como flop
        traza.append([fila[:] + [rhs[i]] for i, fila in enumerate(M)])

    #sustitución regresiva
    
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = rhs[i]
        for j in range(i + 1, n):
            s = s - M[i][j] * x[j]
            cont.retro["mul"] += 1
            cont.retro["add"] += 1
        x[i] = s / M[i][i]
        cont.retro["div"] += 1

    return Factorizacion(M, L, perm, rhs, x, mult_max, crecimiento, traza, pivotes, cont, modo)


# reutilizar L y U

def resolver_con_LU(f, rhs):
    """Resuelve A·x = rhs reutilizando L, U y la permutación (sin refactorizar)."""
    n = len(rhs)
    y = [rhs[f.perm[i]] for i in range(n)]
    for k in range(n - 1):
        for i in range(k + 1, n):
            y[i] = y[i] - f.L[i][k] * y[k]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = y[i]
        for j in range(i + 1, n):
            s = s - f.U[i][j] * x[j]
        x[i] = s / f.U[i][i]
    return x


def resolver_transpuesta_con_LU(f, rhs):
    """Resuelve Aᵀ·x = rhs con la misma factorización (Aᵀ = Uᵀ Lᵀ P)."""
    n = len(rhs)
    z = [0.0] * n
    for i in range(n):                            # Uᵀ z = rhs  (Uᵀ es triangular inferior)
        s = rhs[i]
        for j in range(i):
            s = s - f.U[j][i] * z[j]
        z[i] = s / f.U[i][i]
    w = z[:]
    for i in range(n - 1, -1, -1):                # Lᵀ w = z  (Lᵀ triangular superior, diag 1)
        s = w[i]
        for j in range(i + 1, n):
            s = s - f.L[j][i] * w[j]
        w[i] = s
    x = [0.0] * n
    for i in range(n):                            # x = Pᵀ w
        x[f.perm[i]] = w[i]
    return x
