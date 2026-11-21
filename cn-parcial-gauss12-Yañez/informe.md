# El pivote — Informe técnico

## 1. Identificación

| Campo |
| --- | --- |
| Estudiante | Maximiliano Emmanuel Yáñez Álvarez |
| Sección | 2 |
| Asignatura | Computación Numérica (INF-321), Ingeniería Civil Informática, Universidad Católica del Maule |
| Docentes | Juan Manuel Cabezas Contardo · Pedro Pesante |
| Fecha | Talca, 20/11/2026|


## 2. Planteamiento

Se resuelve A·x = b con n = 12, A densa (144 de 144 elementos no nulos, 9 cifras por elemento) y b con 20 cifras. Los datos se ingresaron como cadenas decimales sin truncar: a doble precisión (`float`, redondeo al valor más cercano, error relativo ≤ 2⁻⁵³ por dato) para el algoritmo, y como racionales exactos (`fractions.Fraction`) para el residuo y la referencia. Truncar b a 7 u 8 cifras no ocurre en ningún punto.

Tres rasgos hacen costoso el problema:

- **Escalas de fila incompatibles.** El máximo absoluto por fila va de 8.05×10⁻⁵ (fila 5) a 9.67×10⁵ (fila 6): diez órdenes de magnitud.

- **Pivote inicial minúsculo.** a(1,1) = 1.3407×10⁻⁸, mientras el mayor elemento de la primera columna vale 3.69×10⁵ (fila 6).

- **Casi dependencia lineal.** Las filas 3, 5 y 7 están muy cerca de ser linealmente dependientes: al ajustar la fila 7 como α·fila 3 + β·fila 5 por mínimos cuadrados (α ≈ 4.177, β ≈ 583.2) el residuo relativo es 2.28×10⁻³ (norma 2; `resultados/residuo.txt`). Es una dependencia del trío {3, 5, 7}: ajustar la fila 5 con las filas 3 y 7 da un residuo parecido (2.0×10⁻³). Esto, más las escalas de fila, da κ₂(A) ≈ 3.08×10¹³.

**Referencia exacta.** Resolver el sistema en aritmética racional da una solución con exactamente 10 decimales (por ejemplo x₁ = 0.7384915627, x₁₂ = −1.0849371625). Es la solución del sistema tal como está dado y se usa solo para medir errores.

## 3. Estrategia de eliminación

Se triangulariza la matriz aumentada [A | b] en 11 pasos (k = 1, …, 11). En cada paso: (1) búsqueda de pivote entre las filas k..12 y eventual intercambio; (2) un multiplicador por fila, m(i,k) = a(i,k)/a(k,k), calculado **una sola vez** y reutilizado en las columnas k+1..12 y en b; (3) actualización de la submatriz solo desde la columna k+1 (el elemento a(i,k) se fija en 0 por construcción, sin gastar flops); (4) al final, sustitución regresiva. Las filas ya fijadas no se vuelven a tocar.

L se guarda (multiplicadores) junto con U y la permutación, de modo que el refinamiento iterativo y el estimador de κ reutilizan la **misma factorización**: no se factoriza dos veces y nunca se forma A⁻¹.

La eliminación densa cuesta 2n³/3 flops en el término dominante (1152 para n = 12) más términos de orden n² (actualizar b y sustituir). El conteo medido (sección 4) es 1354, es decir, no hay trabajo redundante: la variante que recalcula el multiplicador por columna, también instrumentada, cuesta 572 flops más.

## 4. Conteo de operaciones

Convención: 1 flop = 1 operación aritmética elemental (+, −, ×, ÷). `abs()` e intercambios no son flops; las comparaciones de pivoteo y las divisiones del pivoteo escalado se contabilizan aparte. Los números salen de un `Contador` incrementado dentro de los bucles de `src/gauss.py` y se contrastan con la fórmula cerrada solo como control.

| Paso k | Divisiones (multiplicadores) | Multiplicaciones | Sumas/restas | Total flops |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 11 | 132 | 132 | 275 |
| 2 | 10 | 110 | 110 | 230 |
| 3 | 9 | 90 | 90 | 189 |
| 4 | 8 | 72 | 72 | 152 |
| 5 | 7 | 56 | 56 | 119 |
| 6 | 6 | 42 | 42 | 90 |
| 7 | 5 | 30 | 30 | 65 |
| 8 | 4 | 20 | 20 | 44 |
| 9 | 3 | 12 | 12 | 27 |
| 10 | 2 | 6 | 6 | 14 |
| 11 | 1 | 2 | 2 | 5 |
| **Eliminación (11 pasos)** | **66** | **572** | **572** | **1210** |
| Sustitución regresiva | 12 | 66 | 66 | 144 |
| **TOTAL (sin pivoteo)** | **78** | **638** | **638** | **1354** |

Total base: **1354 flops** (716 multiplicaciones/divisiones y 638 sumas/restas). Es el mismo con o sin pivoteo.

**Sobrecosto de estabilidad (aparte).** Pivoteo parcial escalado: 275 operaciones = 132 comparaciones para las escalas sᵢ + 77 divisiones de razón |aᵢₖ|/sᵢ + 66 comparaciones de búsqueda, más 8 intercambios de fila (no son flops). Pivoteo parcial simple: 66 comparaciones y 9 intercambios.

| Estrategia | Flops medidos | Cota de la prueba | Desviación | ¿Dentro de ±10 %? |
| --- | ---: | ---: | ---: | :---: |
| Eliminación densa, multiplicadores reutilizados | 1354 | ≈ 1486 | -8.9 % | sí |
| Lo mismo + pivoteo (búsqueda y escalas) | 1629 | ≈ 1552 (+66 cmp) | +5.0 % | sí |
| Multiplicadores recalculados por columna (variante penalizada) | 1926 | ≈ 1750 | +10.1 % | exceso: +572 flops sobre el mínimo |
| Vía matriz inversa (no implementada, descalifica C2) | – | ≈ 7000 o más | – | – |

Frente a la cota de la estrategia declarada (eliminación con multiplicadores reutilizados + pivoteo, 1486 + 66 = 1552), el costo medido (1354 + 275 = 1629) queda +5.0 %, dentro de ±10 %. La cota de 1486 equivale a calcular explícitamente a(i,k) − m·a(k,k) en la columna que se anula (66 multiplicaciones y 66 restas más); aquí ese cálculo se omite, y por eso el total base queda 8.9 % bajo la referencia.

## 5. Análisis de mantisa y pivoteo

Mismo código, tres corridas (`resultados/residuo.txt`):

| Corrida | Multiplicador máx. | log₂(mult.) | Factor de crecimiento | ‖r‖∞ | Error relativo (norma ∞) | Dígitos mín. correctos | Bits perdidos |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Sin pivoteo | 2.7526×10¹³ | 44.6 | 2.6412×10³ | 1.772×10⁻⁷ | 2.0803×10⁻¹⁰ | 7 | 20.8 |
| Pivoteo parcial | 0.62 | −0.7 | 1.0000 | 5.373×10⁻¹⁰ | 5.5201×10⁻¹³ | 9 | 12.3 |
| Pivoteo parcial escalado | 2.4656×10⁹ | 31.2 | 1.7832 | 2.291×10⁻¹⁰ | 2.3573×10⁻¹³ | 10 | 11.1 |

Factor de crecimiento = máx|a_ij⁽ᵏ⁾| sobre todos los pasos / máx|a_ij⁽⁰⁾|. Bits perdidos = log₂(error relativo / 2⁻⁵³), medidos contra la referencia exacta.

**Por qué a(1,1) es inaceptable.** Sin pivoteo, el pivote es 1.3407×10⁻⁸ y el multiplicador de la fila 6 es 3.69×10⁵ / 1.3407×10⁻⁸ = 2.75×10¹³ ≈ 2⁴⁴·⁶. En la operación a(6,j) − m·a(1,j), el término m·a(1,j) vale entre 7×10⁸ y 2.6×10⁹, mientras a(6,j) está entre 2×10⁵ y 10⁶: el dato original de la fila 6 es solo 2×10⁻⁴ a 10⁻³ del término dominante, de modo que **10 a 12 bits de a(6,j) quedan bajo el último bit** de la mantisa en una sola operación (en la traza del paso 1 sin pivoteo, la fila 6 queda con entradas de orden 10⁹). Ese ruido se propaga durante los 10 pasos siguientes: la matriz activa crece por un factor 2641 (≈ 11.4 bits más) y el error final es 2.1×10⁻¹⁰, es decir 20.8 de los 53 bits perdidos; quedan 32 bits ≈ 9.7 dígitos en norma y 7 dígitos en las componentes pequeñas. El sistema sí se resuelve (ningún pivote se anula), pero con 6 a 7 dígitos en vez de 12.

**Pivoteo parcial frente a escalado.** El parcial simple elige el mayor |a(i,k)| y mantiene todos los multiplicadores ≤ 1; el factor de crecimiento global es 1.00 porque lo domina el elemento 9.67×10⁵ de la fila 6. Aun así, con filas en escalas tan distintas el criterio sin escalar compara magnitudes que no son comparables. El escalado elige el mayor |a(i,k)|/sᵢ (sᵢ = máx|a_ij| original de la fila i). Con esa regla se cumple |a(i,k)|/|a(k,k)| ≤ sᵢ/sₖ, de modo que el término m·a(k,j) restado a la fila i nunca supera la escala sᵢ de esa fila: puede haber multiplicadores enormes (2.5×10⁹) sin crecimiento relativo a la fila (1.78) ni pérdida de bits de la fila destino. Por eso el multiplicador máximo **no** es la medida de estabilidad en este caso; sí lo son el crecimiento y el error medido. El escalado reduce el error en un factor 2.3 respecto al parcial simple y sube el mínimo de 9 a 10 dígitos.

## 6. Condicionamiento

| Magnitud | Valor |
| --- | ---: |
| κ₁(A) estimado (Hager sobre L y U, sin formar A⁻¹) | 4.4675×10¹³ |
| κ₂(A) (control, `numpy.linalg.cond`) | 3.0847×10¹³ |
| ε (doble precisión, 2⁻⁵²) | 2.2204×10⁻¹⁶ |
| Cota κ₂·ε sobre ‖δx‖/‖x‖ | 6.85×10⁻³ (2.2 dígitos) |
| Error real (pivoteo escalado) | 2.36×10⁻¹³ (12.6 dígitos) |
| κ₂ tras equilibrar filas (cada fila / su máx.) | 8.13×10³ |
| κ_eq·ε | 1.81×10⁻¹² (11.7 dígitos) |

La cota clásica ‖δx‖/‖x‖ ≲ κ·ε no garantiza ni tres dígitos, pero el error real es 2.4×10⁻¹³. La cota es pesimista aquí por una razón concreta: casi todo el κ proviene de la escala de las filas, no de la geometría del sistema. Multiplicar una fila por una constante no cambia la solución, y un método con pivoteo escalado es insensible a esa escala; al equilibrar las filas, κ₂ baja de 3.1×10¹³ a 8.1×10³ (el resto, del orden de 10³, corresponde a la casi dependencia de las filas 3, 5 y 7). Con ese número, κ_eq·ε = 1.8×10⁻¹² predice 11.7 dígitos, coherente con los 12.6 observados. Conclusión: en doble precisión son recuperables unos 12 dígitos con pivoteo escalado; la cota κ·ε describe el peor caso, no este.

## 7. Resultados

Solución con pivoteo parcial escalado. «Dígitos afirmados» sale de la primera corrección del refinamiento, |d₁ᵢ|/|xᵢ| (estima el error sin conocer la solución exacta); «verificados» se mide contra la referencia racional. Definición: d dígitos correctos significa error relativo < 10⁻ᵈ (el valor impreso queda a menos de una unidad de su última cifra). Coinciden en las 12 incógnitas.

| i | x reportada (sin refinar) | Dígitos afirmados | Dígitos verificados vs. exacta | x refinada (17 cifras) |
| ---: | ---: | ---: | ---: | ---: |
| x1 | 7.384915626999e-01 | 13 | 13 | 7.3849156270000005e-01 |
| x2 | -1.20593741680e+00 | 12 | 12 | -1.2059374168000001e+00 |
| x3 | 4.387195200e-03 | 10 | 10 | 4.3871952000000001e-03 |
| x4 | -9.1726453800e-02 | 11 | 11 | -9.1726453799999996e-02 |
| x5 | 2.84610937250e+00 | 12 | 12 | 2.8461093725000000e+00 |
| x6 | -6.281540899e-03 | 10 | 10 | -6.2815409000000003e-03 |
| x7 | 1.47382059160e+00 | 12 | 12 | 1.4738205916000000e+00 |
| x8 | -3.8492716400e-02 | 11 | 11 | -3.8492716400000002e-02 |
| x9 | 9.15264803700e-01 | 12 | 12 | 9.1526480369999996e-01 |
| x10 | -2.60713948520e+00 | 12 | 12 | -2.6071394851999998e+00 |
| x11 | 5.1873629400e-02 | 11 | 11 | 5.1873629400000003e-02 |
| x12 | -1.08493716250e+00 | 12 | 12 | -1.0849371624999999e+00 |

El mínimo es 10 dígitos (x₃ y x₆, componentes pequeñas); la norma da 12.6.

## 8. Verificación

**Residuo en precisión extendida.** r = b − A·x se calcula con A, b racionales exactos y x tal cual (sin redondeo). Con el pivoteo escalado, ‖r‖∞ = 2.291×10⁻¹⁰ (‖b‖∞ = 7.81×10⁵). El residuo de cada fila es proporcional a su escala (entre 10⁻²⁰ y 10⁻¹⁷ en las filas pequeñas, entre 10⁻¹² y 10⁻¹⁰ en las grandes). El error inverso relativo ‖r‖∞/(‖A‖∞‖x‖∞) = 1.1×10⁻¹⁷, menor que la unidad de redondeo u = 2⁻⁵³ ≈ 1.1×10⁻¹⁶: el algoritmo es estable hacia atrás y el error restante viene del condicionamiento.

**Refinamiento iterativo** (reutiliza L y U, residuo exacto, corrección en doble; costo ≈ 276 flops por pasada):

| Pasada | ‖r‖∞ antes | ‖d‖∞/‖x‖∞ | Error relativo, iterado en doble | Error relativo, iterado extendido |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 2.291×10⁻¹⁰ | – | 2.36×10⁻¹³ | 2.36×10⁻¹³ |
| 1 | 2.291×10⁻¹⁰ | 2.36×10⁻¹³ | 6.02×10⁻¹⁷ | 6.17×10⁻²⁶ |
| 2 | 2.603×10⁻¹⁰ | 6.02×10⁻¹⁷ | 6.02×10⁻¹⁷ | 1.61×10⁻³⁸ |
| 3 | 2.603×10⁻¹⁰ | 6.02×10⁻¹⁷ | 6.02×10⁻¹⁷ | 4.22×10⁻⁵¹ |

Con el iterado guardado en doble, una pasada lleva el error a 6×10⁻¹⁷ (la mitad de una unidad en el último lugar de x) y ahí se detiene: el residuo de un x en doble no baja de ≈ u·‖A‖·‖x‖ por el redondeo del propio x. Con el iterado acumulado en precisión extendida la convergencia es lineal, con un factor ≈ 10⁻¹³ por pasada, lo que confirma que L y U sirven como inversa aproximada pese a κ ≈ 3×10¹³.

**Control con biblioteca** (no sustituye la implementación propia): `numpy.linalg.solve` da error relativo 5.2×10⁻¹³ (9 dígitos mínimos); la diferencia máxima con la solución propia es 2.9×10⁻¹³ (relativa a ‖x‖∞), del mismo orden de ruido que el error de ambas.

## 9. Conclusiones

El compromiso elegido es **eliminación densa sin trabajo redundante (1354 flops) más pivoteo parcial escalado (+275 operaciones, +20 %)**. Ese sobrecosto es de orden n² y no cambia el costo cúbico; a cambio, el error de norma baja de 2.1×10⁻¹⁰ a 2.4×10⁻¹³ (unas 880 veces) y el mínimo sube de 7 a 10 dígitos. Una pasada de refinamiento (≈ 276 flops, sin refactorizar) lleva el resultado al límite de la doble precisión. El camino más barato (sin pivoteo) no es el más estable, y el más estable (pivoteo + refinamiento) cuesta apenas un 20–40 % más: aquí el conflicto del enunciado se resuelve a favor de la estabilidad.

Con n = 1000: la eliminación cuesta 2n³/3 ≈ 6.7×10⁸ flops, y la búsqueda de pivote y las escalas (≈ 1.5×10⁶ operaciones) pasan a ser despreciables. Usaría LU por bloques con pivoteo parcial (nivel 3 de BLAS, como LAPACK `dgetrf`) por rendimiento de caché, equilibraría las filas una vez antes de factorizar, estimaría κ con un estimador tipo `dgecon` y haría refinamiento iterativo con residuo en precisión extendida o con suma compensada. Nunca invertiría A.

## 10. Reproducibilidad

- Python 3.12.3; numpy 2.5.4; mpmath 1.3.0 (`requirements.txt`). El algoritmo propio solo usa la biblioteca estándar (`fractions`); numpy se usa para el control (`linalg.solve`, `linalg.cond`) y mpmath para imprimir cifras extendidas.
- Comando: `pip install -r requirements.txt` y luego `python src/main.py`; regenera todo `resultados/` sin intervención manual.

**Fuentes.**

1. N. J. Higham, *Accuracy and Stability of Numerical Algorithms*, 2.ª ed., SIAM, 2002 (pivoteo, factor de crecimiento, error inverso).
2. G. H. Golub y C. F. Van Loan, *Matrix Computations*, 4.ª ed., Johns Hopkins University Press, 2013 (costo de la eliminación, equilibrado de filas).
3. W. W. Hager, «Condition estimates», *SIAM J. Sci. Stat. Comput.* 5(2), 311–316, 1984 (estimador de ‖A⁻¹‖₁).
4. IEEE Std 754-2019, *IEEE Standard for Floating-Point Arithmetic*.
