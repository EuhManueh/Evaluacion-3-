# cn-parcial-gauss12

Resuelve un sistema lineal denso 12×12 por eliminación gaussiana propia (con y sin pivoteo), con conteo de flops y verificación en precisión extendida.

## Autoría y asignatura

Maximiloiano Yañez Alvarez — Computación Numérica (INF-321), Ingeniería Civil Informática, Universidad Católica del Maule. Prueba Parcial «El pivote», 20/11/2026
## Requisitos

Python 3.10 o superior (probado con 3.12.3).

```
pip install -r requirements.txt
```

## Ejecución

```
python src/main.py
```

Un solo comando, sin pasos manuales. Escribe en `resultados/`: `solucion.csv`, `traza_eliminacion.txt`, `conteo_flops.md` y `residuo.txt`. También imprime un resumen por pantalla.

## Estructura de carpetas

```
datos/        matriz A (12x12) y vector b con todas sus cifras
src/          gauss.py (eliminación y contador), conteo.py, verificacion.py, main.py
resultados/   salidas generadas por main.py
informe.md    informe técnico y análisis de error
```

## Resultado resumido

Incógnitas a 10 dígitos significativos (pivoteo parcial escalado y una pasada de refinamiento):

| Incógnita | Valor (10 dígitos) | Dígitos afirmados (sin refinar) |
| --- | ---: | ---: |
| x1 | 7.384915627e-01 | 13 |
| x2 | -1.205937417e+00 | 12 |
| x3 | 4.387195200e-03 | 10 |
| x4 | -9.172645380e-02 | 11 |
| x5 | 2.846109372e+00 | 12 |
| x6 | -6.281540900e-03 | 10 |
| x7 | 1.473820592e+00 | 12 |
| x8 | -3.849271640e-02 | 11 |
| x9 | 9.152648037e-01 | 12 |
| x10 | -2.607139485e+00 | 12 |
| x11 | 5.187362940e-02 | 11 |
| x12 | -1.084937162e+00 | 12 |

- Flops de eliminación y sustitución: **1354** (+ 275 de sobrecosto del pivoteo escalado).
- Norma infinito del residuo exacto, sin refinar: **‖r‖∞ = 2.291×10⁻¹⁰**.

## Más detalle

El análisis completo está en [informe.md](informe.md).

## Licencia / nota de uso

Trabajo académico de la UCM. Se permite su lectura y ejecución con fines de evaluación.
