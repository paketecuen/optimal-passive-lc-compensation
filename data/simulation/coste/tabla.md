# Coste de calculo en el mismo host

Generado desde [comparacion.csv](comparacion.csv). Hipotesis, metodos y limites: [guia](../../docs/COMPARACION_PARA_EL_PAPER.md). Los guiones representan datos no disponibles, no ceros.

Medidas MATLAB `24.1.0.2568132 (R2024a) Update 1`, CPU `Intel(R) Core(TM) i9-9880H CPU @ 2.30GHz`, mismo registro y equipo. `timeit` estima el tiempo mediante llamadas repetidas; no es la mediana de tres diseños de la tabla estacionaria. Tiempos orientativos; memoria de arrays estimada, sin runtime ni temporales internos.

## Cuentas de operaciones y memoria

$N$: muestras, $M$: ramas, $S$: estados admisibles, $H$: coeficientes del soporte conocido. Para un banco binario, $S=2^M$. Las plantillas, su Gram y las energías por estado se preparan antes de medir.

| Operación | Trabajo dominante | Memoria principal |
|---|---|---|
| Preparar $G=\Phi^T\Phi/N$ | $O(NM^2)$ | $NM+M^2$ reales |
| Correlaciones y energía por ventana | $M(2N-1)+2N-1+O(M)$ operaciones reales | $NM+N+M$ reales |
| Evaluar $S$ costes con energías preparadas | Aproximadamente $2MS+O(S)$ | Tabla de estados $MS$ y energías $S$ |
| FFT completa y correlaciones | $O(N\log N+NM+MS)$ | $NM+N$ complejos y $MS+S$ reales |
| FFT con soporte conocido | $O(N\log N+HM+MS)$ | $HM+N$ complejos y $MS+S$ reales |
| Enumeración directa vectorizada | $O(NMS)$ | Plantillas $NM$, estados $MS$, intermedio $NS$ |

La cuenta de FFT depende de la factorización de N y la biblioteca; no se le atribuye una constante universal. No se incluyen diseño/adaptación de plantillas, estimación de frecuencia, filtrado ni ejecutor. La dimensión del espacio de estados sigue siendo exponencial. Estas medidas no prueban que cualquier implementación temporal sea más rápida que una FFT optimizada.

## Medidas en MATLAB

| N | M | Estados | Metodo | Tiempo_s | Memoria_estimacion_bytes | Error_coste_A2 |
|---|---|---|---|---|---|---|
| 63 | 2 | 4 | Temporal plantillas | 1.1772755e-06 | 1640 | 0 |
| 63 | 2 | 4 | FFT completa | 6.249306e-06 | 3120 | 1.7763568e-15 |
| 63 | 2 | 4 | FFT soporte conocido | 4.6461097e-06 | 1296 | 1.7763568e-15 |
| 63 | 2 | 4 | Enumeracion directa | 1.7963371e-06 | 3592 | 2.8421709e-14 |
| 63 | 4 | 16 | Temporal plantillas | 1.105735e-06 | 3288 | 0 |
| 63 | 4 | 16 | FFT completa | 5.7100989e-06 | 5680 | 7.1054274e-15 |
| 63 | 4 | 16 | FFT soporte conocido | 3.9029481e-06 | 2032 | 7.1054274e-15 |
| 63 | 4 | 16 | Enumeracion directa | 2.8628163e-06 | 11096 | 8.5265128e-14 |
| 63 | 8 | 256 | Temporal plantillas | 1.4783667e-06 | 23480 | 0 |
| 63 | 8 | 256 | FFT completa | 5.7087008e-06 | 27504 | 5.6843419e-14 |
| 63 | 8 | 256 | FFT soporte conocido | 4.9456755e-06 | 20208 | 5.6843419e-14 |
| 63 | 8 | 256 | Enumeracion directa | 5.1444045e-05 | 149944 | 1.1368684e-13 |
| 255 | 2 | 4 | Temporal plantillas | 1.0039178e-06 | 6248 | 0 |
| 255 | 2 | 4 | FFT completa | 8.5144278e-06 | 12336 | 1.7763568e-15 |
| 255 | 2 | 4 | FFT soporte conocido | 6.9149722e-06 | 4368 | 1.7763568e-15 |
| 255 | 2 | 4 | Enumeracion directa | 2.5304784e-06 | 14344 | 1.4210855e-14 |
| 255 | 4 | 16 | Temporal plantillas | 1.228773e-06 | 10968 | 0 |
| 255 | 4 | 16 | FFT completa | 1.0206316e-05 | 21040 | 1.7763568e-15 |
| 255 | 4 | 16 | FFT soporte conocido | 7.6054477e-06 | 5104 | 1.7763568e-15 |
| 255 | 4 | 16 | Enumeracion directa | 5.2535636e-06 | 43352 | 8.5265128e-14 |
| 255 | 8 | 256 | Temporal plantillas | 2.0849677e-06 | 37304 | 0 |
| 255 | 8 | 256 | FFT completa | 1.155765e-05 | 55152 | 5.6843419e-14 |
| 255 | 8 | 256 | FFT soporte conocido | 8.7281514e-06 | 23280 | 5.6843419e-14 |
| 255 | 8 | 256 | Enumeracion directa | 6.9907986e-05 | 556984 | 1.4210855e-13 |
| 1023 | 2 | 4 | Temporal plantillas | 2.671809e-06 | 24680 | 0 |
| 1023 | 2 | 4 | FFT completa | 7.1848016e-05 | 49200 | 1.4210855e-14 |
| 1023 | 2 | 4 | FFT soporte conocido | 3.8483316e-05 | 16656 | 1.4210855e-14 |
| 1023 | 2 | 4 | Enumeracion directa | 1.0615439e-05 | 57352 | 7.1054274e-14 |
| 1023 | 4 | 16 | Temporal plantillas | 4.598496e-06 | 41688 | 0 |
| 1023 | 4 | 16 | FFT completa | 7.614436e-05 | 82480 | 2.8421709e-14 |
| 1023 | 4 | 16 | FFT soporte conocido | 4.1245339e-05 | 17392 | 2.8421709e-14 |
| 1023 | 4 | 16 | Enumeracion directa | 3.2139476e-05 | 172376 | 1.1368684e-13 |
| 1023 | 8 | 256 | Temporal plantillas | 4.7821371e-06 | 92600 | 0 |
| 1023 | 8 | 256 | FFT completa | 5.6794104e-05 | 165744 | 5.6843419e-14 |
| 1023 | 8 | 256 | FFT soporte conocido | 3.234047e-05 | 35568 | 5.6843419e-14 |
| 1023 | 8 | 256 | Enumeracion directa | 0.00020576843 | 2185144 | 1.7053026e-13 |
