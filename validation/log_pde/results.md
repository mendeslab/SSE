# Measured log-PDE results

Successful runs: 183/183.

Exact exterior values; fixed rtol=1e-12; single-run Python timings.

## Spatial errors at atol=1e-9 and domain [−4,4]

| Case | Order | Δx=0.5 | 0.25 | 0.125 | 0.0625 | Final observed order |
|---|---:|---:|---:|---:|---:|---:|
| narrow-centered | 2 | 1.6e-12 | 1.76e-13 | 4.4e-15 | 5.55e-16 | quadratic log; no spatial-order inference |
| narrow-centered | 4 | 1.03e-12 | 2.5e-13 | 6.46e-15 | 8.08e-16 | quadratic log; no spatial-order inference |
| narrow-centered | 6 | 3.5e-12 | 2.3e-12 | 2.1e-14 | 3.13e-16 | quadratic log; no spatial-order inference |
| narrow-offset | 2 | 1.38e-12 | 2e-13 | 4.86e-15 | 7.35e-16 | quadratic log; no spatial-order inference |
| narrow-offset | 4 | 4.46e-13 | 2.9e-13 | 9.12e-15 | 5.69e-16 | quadratic log; no spatial-order inference |
| narrow-offset | 6 | 2.4e-12 | 2.53e-12 | 2.76e-14 | 1.15e-15 | quadratic log; no spatial-order inference |
| broad | 2 | 3.88e-13 | 7.75e-14 | 1.03e-14 | 3.81e-15 | quadratic log; no spatial-order inference |
| broad | 4 | 4.05e-13 | 3.14e-14 | 7.47e-15 | 3.83e-15 | quadratic log; no spatial-order inference |
| broad | 6 | 6.5e-13 | 2.95e-14 | 6.98e-15 | 4.47e-15 | quadratic log; no spatial-order inference |
| overlapping | 2 | 0.0209 | 0.00617 | 0.00151 | 0.000376 | 2.01 |
| overlapping | 4 | 0.00939 | 0.000833 | 5.58e-05 | 3.63e-06 | 3.94 |
| overlapping | 6 | 0.0118 | 0.000405 | 7.96e-06 | 1.4e-07 | 5.83 |
| separated | 2 | 0.0391 | 0.0203 | 0.00521 | 0.0013 | 2.00 |
| separated | 4 | 0.142 | 0.0179 | 0.00125 | 7.87e-05 | 3.99 |
| separated | 6 | 0.214 | 0.0215 | 0.000549 | 1.06e-05 | 5.70 |
| merge | 2 | 0.166 | 0.0399 | 0.0104 | 0.00262 | 1.99 |
| merge | 4 | 0.115 | 0.011 | 0.00085 | 5.65e-05 | 3.91 |
| merge | 6 | 0.1 | 0.00532 | 0.000153 | 2.93e-06 | 5.71 |

## Domain and tight-tolerance checks

| Case | Extent | Duration | atol | Grid L2 | Interpolated L2 | Integral error | Quadrature change 16→32 | Evaluations | Seconds |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| narrow-centered | 4 | 0.2 | 1e-06 | 3.61e-11 | 5.46e-11 | 2.5e-12 | 4.44e-16 | 986 | 0.060 |
| narrow-centered | 4 | 0.2 | 1e-09 | 2.1e-14 | 2.94e-14 | 3.55e-15 | 4.44e-16 | 3530 | 0.213 |
| narrow-offset | 4 | 0.2 | 1e-06 | 4.66e-11 | 6.45e-11 | 9.26e-13 | 4.44e-16 | 992 | 0.061 |
| narrow-offset | 4 | 0.2 | 1e-09 | 2.76e-14 | 3.54e-14 | 0 | 4.44e-16 | 3548 | 0.215 |
| broad | 4 | 0.2 | 1e-06 | 5.19e-12 | 4.61e-12 | 9.7e-13 | 0 | 26 | 0.002 |
| broad | 4 | 0.2 | 1e-09 | 6.98e-15 | 6.22e-15 | 8.16e-16 | 0 | 74 | 0.005 |
| overlapping | 4 | 0.2 | 1e-06 | 7.96e-06 | 7.95e-06 | 3.02e-06 | 0 | 62 | 0.004 |
| overlapping | 4 | 0.2 | 1e-09 | 7.96e-06 | 7.95e-06 | 3.02e-06 | 0 | 212 | 0.014 |
| separated | 4 | 0.2 | 1e-06 | 0.000549 | 0.000649 | 6.38e-05 | 0 | 128 | 0.008 |
| separated | 4 | 0.2 | 1e-09 | 0.000549 | 0.000649 | 6.38e-05 | 0 | 470 | 0.030 |
| merge | 4 | 0.2 | 1e-06 | 0.000153 | 0.000153 | 5.64e-05 | 0 | 324 | 0.022 |
| merge | 4 | 0.2 | 1e-09 | 0.000153 | 0.000153 | 5.64e-05 | 0 | 1188 | 0.076 |
| narrow-offset | 2 | 0.2 | 1e-09 | 8.73e-14 | 9.95e-14 | 5.42e-14 | 4.44e-16 | 2018 | 0.122 |
| narrow-offset | 8 | 0.2 | 1e-09 | 2.12e-14 | 2.7e-14 | 0 | 0 | 5804 | 0.357 |
| separated | 2 | 0.2 | 1e-09 | 0.000549 | 0.000649 | 6.39e-05 | 0 | 194 | 0.013 |
| separated | 8 | 0.2 | 1e-09 | 0.000549 | 0.000649 | 6.38e-05 | 0 | 1238 | 0.079 |
| narrow-centered | 4 | 0.2 | 1e-11 | 8.59e-15 | 1.05e-14 | 1.78e-15 | 0 | 4466 | 0.271 |
| narrow-offset | 4 | 0.2 | 1e-11 | 1.09e-14 | 1.25e-14 | 8.88e-16 | 0 | 4478 | 0.274 |
| broad | 4 | 0.2 | 1e-11 | 2.07e-16 | 2.62e-16 | 7.23e-17 | 0 | 146 | 0.009 |
| overlapping | 4 | 0.2 | 1e-11 | 7.96e-06 | 7.95e-06 | 3.02e-06 | 0 | 386 | 0.025 |
| separated | 4 | 0.2 | 1e-11 | 0.000549 | 0.000649 | 6.38e-05 | 0 | 800 | 0.051 |
| merge | 4 | 0.2 | 1e-11 | 0.000153 | 0.000153 | 5.64e-05 | 0 | 2046 | 0.130 |
| width-0.16 | 4 | 0.2 | 1e-09 | 1.44e-15 | 1.92e-15 | 0 | 4.44e-16 | 1466 | 0.089 |
| width-0.32 | 4 | 0.2 | 1e-09 | 2.31e-15 | 2.3e-15 | 8.88e-16 | 0 | 410 | 0.025 |
| narrow-offset | 4 | 0.02 | 1e-06 | 1.29e-13 | 7.33e-12 | 1.11e-12 | 0 | 290 | 0.018 |
| narrow-offset | 4 | 0.02 | 1e-09 | 5.84e-16 | 2.19e-14 | 3.11e-15 | 0 | 968 | 0.061 |
| narrow-offset | 4 | 1 | 1e-06 | 1.21e-11 | 1.21e-11 | 6.99e-12 | 0 | 1496 | 0.091 |
| narrow-offset | 4 | 1 | 1e-09 | 1.04e-14 | 1.04e-14 | 5.33e-15 | 0 | 5450 | 0.334 |
| separated | 4 | 0.02 | 1e-06 | 0.000308 | 0.000726 | 1.23e-05 | 8.88e-16 | 20 | 0.002 |
| separated | 4 | 0.02 | 1e-09 | 0.000308 | 0.000726 | 1.23e-05 | 0 | 68 | 0.004 |
| separated | 4 | 1 | 1e-06 | 0.000151 | 0.000151 | 8.1e-05 | 0 | 326 | 0.021 |
| separated | 4 | 1 | 1e-09 | 0.000151 | 0.000151 | 8.1e-05 | 0 | 1214 | 0.077 |

## Matched-accuracy comparisons

Smallest observed propagation time among the measured grids at atol=1e-9, domain [−4,4],
duration 0.2. Target is relative L2 error of the exponentiated log interpolant.
Timings are single-run measurements, including analytic exterior evaluations; no native comparison.

| Case | Target | Order | Δx | Points | Error | Evaluations | Seconds |
|---|---:|---:|---:|---:|---:|---:|---:|
| overlapping | 0.001 | 2 | 0.0625 | 129 | 0.000375 | 224 | 0.014 |
| overlapping | 0.001 | 4 | 0.25 | 33 | 0.000811 | 152 | 0.010 |
| overlapping | 0.001 | 6 | 0.25 | 33 | 0.000391 | 164 | 0.011 |
| overlapping | 0.0001 | 2 | 0.03125 | 257 | 9.38e-05 | 548 | 0.034 |
| overlapping | 0.0001 | 4 | 0.125 | 65 | 5.56e-05 | 194 | 0.012 |
| overlapping | 0.0001 | 6 | 0.125 | 65 | 7.95e-06 | 212 | 0.014 |
| overlapping | 1e-05 | 2 | 0.0078125 | 1025 | 5.86e-06 | 3350 | 0.221 |
| overlapping | 1e-05 | 4 | 0.0625 | 129 | 3.63e-06 | 278 | 0.018 |
| overlapping | 1e-05 | 6 | 0.125 | 65 | 7.95e-06 | 212 | 0.014 |
| separated | 0.001 | 2 | 0.03125 | 257 | 0.000325 | 674 | 0.042 |
| separated | 0.001 | 4 | 0.0625 | 129 | 7.87e-05 | 590 | 0.037 |
| separated | 0.001 | 6 | 0.125 | 65 | 0.000649 | 470 | 0.030 |
| separated | 0.0001 | 2 | 0.015625 | 513 | 8.14e-05 | 1790 | 0.112 |
| separated | 0.0001 | 4 | 0.0625 | 129 | 7.87e-05 | 590 | 0.037 |
| separated | 0.0001 | 6 | 0.0625 | 129 | 1.06e-05 | 638 | 0.042 |
| separated | 1e-05 | 2 | — | — | not reached | — | — |
| separated | 1e-05 | 4 | 0.03125 | 257 | 5.06e-06 | 842 | 0.054 |
| separated | 1e-05 | 6 | 0.03125 | 257 | 1.89e-07 | 944 | 0.063 |
| merge | 0.001 | 2 | 0.03125 | 257 | 0.000656 | 1746 | 0.107 |
| merge | 0.001 | 4 | 0.125 | 65 | 0.00085 | 1116 | 0.069 |
| merge | 0.001 | 6 | 0.125 | 65 | 0.000153 | 1188 | 0.076 |
| merge | 0.0001 | 2 | 0.0078125 | 1025 | 4.1e-05 | 11892 | 0.773 |
| merge | 0.0001 | 4 | 0.0625 | 129 | 5.65e-05 | 1386 | 0.086 |
| merge | 0.0001 | 6 | 0.0625 | 129 | 2.93e-06 | 1482 | 0.095 |
| merge | 1e-05 | 2 | — | — | not reached | — | — |
| merge | 1e-05 | 4 | 0.03125 | 257 | 3.6e-06 | 2088 | 0.134 |
| merge | 1e-05 | 6 | 0.0625 | 129 | 2.93e-06 | 1482 | 0.095 |
