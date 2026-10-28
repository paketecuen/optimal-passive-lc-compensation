# Mismas condiciones: criterios de maniobra

Generado desde [comparacion_condensada.csv](comparacion_condensada.csv). Hipotesis, metodos y limites: [guia](../../docs/COMPARACION_PARA_EL_PAPER.md). Los guiones representan datos no disponibles, no ceros.

Los casos históricos 1, 2 y 3 comparten condiciones físicas y resultados: se presentan una sola vez por criterio. `Casos_historicos` conserva la correspondencia. Las [15 filas originales](comparacion.csv) permanecen intactas. Esta tabla contiene nueve combinaciones de tres condiciones físicas y tres criterios.

| Casos_historicos | Criterio | R_ohm | Rb_ohm | VC0_V | I_pico_A | VC_pico_V | Espera_max_ms | Espera_media_ms | Pendiente_max_A_s | Desajuste_max_mJ | Minimo_valvula_max_V | Residuo_energia_J | Residuo_OFF_A |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1/2/3 | inmediato | 0.05 | 1000 | 0 | 3.6899529 | 5.31687 | 0 | 0 | 1080.1602 | 3.3146878 | 0.00067576877 | 1.5632634e-13 | 9.4312058e-13 |
| 1/2/3 | valvula | 0.05 | 1000 | 0 | 3.9314204 | 5.7938025 | 1.5 | 0.84734051 | 0.00024442564 | 6.3018656 | 0.0022836334 | 1.6946861e-13 | 9.3747059e-13 |
| 1/2/3 | estado | 0.05 | 1000 | 0 | 2.9268116 | 4.1583253 | 18.850801 | 10.635816 | 1475.0265 | 0.46537343 | 0.0027181932 | 1.3994361e-13 | 4.3720106e-12 |
| 4 | inmediato | 0 | ∞ | 0 | 3.9630264 | 5.7335233 | 0 | 0 | 1080.1602 | 2.9718581 | 0.0026786709 | 1.687539e-13 | 6.9967851e-12 |
| 4 | valvula | 0 | ∞ | 0 | 4.4819532 | 6.4335994 | 18.156354 | 9.8281772 | 0.00029115562 | 6.3466158 | 0.0055204774 | 1.6651535e-13 | 1.5574261e-13 |
| 4 | estado | 0 | ∞ | 0 | 3.0877739 | 4.5333159 | 19.000011 | 10.720589 | 1560.6949 | 0.47758143 | 0.00022781437 | 1.5460029e-13 | 1.3256722e-12 |
| 5 | inmediato | 0.05 | 1000 | 6.363961 | 5.2225529 | 7.5704151 | 0 | 0 | 962.90379 | 22.391218 | 1.9629241 | 2.3990532e-13 | 1.9976798e-12 |
| 5 | valvula | 0.05 | 1000 | 6.363961 | 5.080455 | 7.2415838 | 19.002001 | 10.645773 | 638.56316 | 23.455699 | 1.9629241 | 1.9560048e-13 | 1.330186e-14 |
| 5 | estado | 0.05 | 1000 | 6.363961 | 5.0378773 | 7.1296122 | 2.4219764 | 2.4214067 | 2531.4114 | 9.891254 | 1.9629241 | 1.7862795e-13 | 6.4045991e-15 |
