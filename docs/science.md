# NMR Methodology

> This guide will be completed in Phase 30.

## Topics to Be Covered

This guide will explain the scientific methodology behind NMR chemical shift predictions:

- **NMR Chemical Shift Prediction Fundamentals**
  - What chemical shifts measure
  - Why quantum mechanical calculations are needed
  - Comparison with empirical methods

- **DFT Theory**
  - B3LYP functional selection
  - Basis set considerations (6-31G*, 6-311+G(2d,p))
  - GIAO method for magnetic properties
  - Computational cost vs accuracy tradeoffs

- **COSMO Solvation Model**
  - Implicit solvation approach
  - Dielectric constant effects
  - Solvent-dependent shift corrections

- **Linear Scaling Methodology**
  - Reference compound calibration
  - Slope and intercept determination
  - Per-nucleus scaling factors

- **DELTA50 Benchmark**
  - Benchmark dataset description
  - Scaling factor derivation
  - Expected accuracy metrics
  - Comparison with literature methods

- **Boltzmann Weighting**
  - Conformer energy distribution
  - Population-weighted averaging
  - Temperature dependence
  - Ensemble size considerations

- **Expected Accuracy and Limitations**
  - Typical error ranges for 1H and 13C
  - Factors affecting accuracy
  - When predictions may fail

- **Literature References**
  - ISiCLE methodology (Pacific Northwest National Laboratory)
  - DELTA50 benchmark (Grimblat et al.)
  - CREST conformer generation (Grimme group)
  - DP4+ probability analysis (Goodman group)

---

For the original ISiCLE publication, see [PNNL ISiCLE](https://github.com/pnnl/isicle).
