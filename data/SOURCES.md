# Benchmark provenance

Carlos Henrique Marchi, Roberta Suero and Luciano Kiyoshi Araki (2009),
**The Lid-Driven Square Cavity Flow: Numerical Solution with a 1024 x 1024 Grid**,
Journal of the Brazilian Society of Mechanical Sciences and Engineering,
31(3), 186–198. DOI: **10.1590/S1678-58782009000300004**.

- Publisher full text: https://www.scielo.br/j/jbsmse/a/XLCc5LvzmPJTtLCmYYXbkKb/?lang=en
- Published-paper PDF mirror used for visual verification:
  https://scispace.com/pdf/the-lid-driven-square-cavity-flow-numerical-solution-with-a-30hd5rlefk.pdf
- Data location: **Table 6, printed page 194 (PDF page 9)**.
- Checked against the rendered table on 11 October 2026.

`marchi_2009_table6_centerlines.csv` transcribes only the numerical centreline
velocity data for Re=10, 100, 400 and 1000: 15 interior positions per component.
For u, `coordinate` is y at x=0.5; for v it is x at y=0.5. Velocities are
normalised by lid speed. The reference problem is the classical steady,
2D incompressible unit square: uniform unit-speed top lid and stationary
other walls. This is not the paper's separate manufactured-solution problem.

These Table 6 velocity values are extrapolated results based on the
1024 x 1024 study, not raw coarse-grid values and not experimental data.
The article provides estimated discretisation errors in Table 7. They are
reference numerical values, not exact analytical solutions. We compare our
interpolated staggered-grid profiles at exactly the 15 tabulated coordinates.
No artificial wall endpoints are included in RMSE or maximum-error metrics.

The paper's mesh counts/discretisation and coupling are not identical to ours:
this is a same-physical-case validation, not a same-method timing comparison.
Only numerical data are redistributed here; the article PDF is not bundled.
