# Conteo de operaciones (n = 12)

**Convención.** 1 flop = 1 operación aritmética elemental (+, −, ×, ÷) en punto flotante. Los números salen de `Contador`, que se incrementa dentro de los bucles de `src/gauss.py`. El elemento a(i,k) que se anula no se calcula (se fija en 0); cada multiplicador se calcula una sola vez y se reutiliza en toda la fila y en b. `abs()` y los intercambios de fila no se cuentan como flops (los intercambios se informan aparte). Si además se calculase explícitamente a(i,k) − m·a(k,k) (66 multiplicaciones y 66 restas más) el total sería 1486, que es la cota de referencia de la prueba.

## Eliminación + sustitución regresiva (costo base, idéntico con o sin pivoteo)

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

Control cruzado: la fórmula cerrada da 716 multiplicaciones/divisiones y 638 sumas/restas (1354 en total) y coincide con el contador.

## Sobrecosto de estabilidad: pivoteo parcial escalado (contado aparte)

| Paso k | Divisiones de razón (escalado) | Comparaciones | Intercambios de fila |
| ---: | ---: | ---: | ---: |
| 1 | 12 | 11 | 1 |
| 2 | 11 | 10 | 0 |
| 3 | 10 | 9 | 1 |
| 4 | 9 | 8 | 1 |
| 5 | 8 | 7 | 1 |
| 6 | 7 | 6 | 1 |
| 7 | 6 | 5 | 1 |
| 8 | 5 | 4 | 1 |
| 9 | 4 | 3 | 1 |
| 10 | 3 | 2 | 0 |
| 11 | 2 | 1 | 0 |
| Escalas s_i (máx por fila) | – | 132 | – |
| **Sobrecosto total de estabilidad** | **77** | **198** | **8** |

Pivoteo parcial simple (sin escalar): 66 comparaciones, 9 intercambios, sin divisiones ni escalas.

## Comparación con las cotas de la prueba

| Estrategia | Flops medidos | Cota de la prueba | Desviación | ¿Dentro de ±10 %? |
| --- | ---: | ---: | ---: | :---: |
| Eliminación densa, multiplicadores reutilizados | 1354 | ≈ 1486 | -8.9 % | sí |
| Lo mismo + pivoteo (búsqueda y escalas) | 1629 | ≈ 1552 (+66 cmp) | +5.0 % | sí |
| Multiplicadores recalculados por columna (variante penalizada) | 1926 | ≈ 1750 | +10.1 % | exceso: +572 flops sobre el mínimo |
| Vía matriz inversa (no implementada, descalifica C2) | – | ≈ 7000 o más | – | – |
