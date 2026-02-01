# NMR Methodology

This guide explains the scientific methodology underlying NMR chemical shift predictions in the QM NMR Calculator. It provides the theoretical foundation for both academic researchers seeking derivations and literature references, and developers needing to understand the implementation.

---

## Introduction

### What is NMR Chemical Shift Prediction?

Nuclear Magnetic Resonance (NMR) spectroscopy measures the magnetic environment of atomic nuclei in a molecule. The **chemical shift** observed for each nucleus depends on its electronic environment, making NMR a powerful tool for structure elucidation. However, interpreting experimental spectra---especially for novel compounds or complex stereochemistry---can be challenging.

**Computational NMR chemical shift prediction** uses quantum mechanical calculations to predict what chemical shifts a proposed structure *should* exhibit. By comparing predicted shifts to experimental data, researchers can:

- **Verify proposed structures:** Does the predicted spectrum match what was observed?
- **Assign stereochemistry:** Which diastereomer's predictions better match experiment?
- **Guide synthesis:** What shifts should we expect for a target molecule?
- **Resolve ambiguities:** When multiple structures could fit the data, predictions can discriminate

### Why Computational Prediction?

Empirical methods (increment rules, group additivity) fail for complex molecules because chemical shifts depend on subtle electronic effects that propagate through molecular structure. Only quantum mechanical calculations can accurately capture:

- Through-bond electronic effects (inductive, mesomeric)
- Through-space interactions (anisotropic shielding)
- Conformational averaging
- Solvent effects on electron density

### Methodology Pipeline

The QM NMR Calculator implements a multi-stage pipeline for chemical shift prediction:

```
                    +------------------+
                    | Input Structure  |
                    |    (SMILES)      |
                    +--------+---------+
                             |
                             v
                 +-----------------------+
                 | Conformer Generation  |
                 |  (RDKit or CREST)     |
                 +-----------+-----------+
                             |
                             v
                 +-----------------------+
                 | Geometry Optimization |
                 |  (B3LYP/6-31G*)       |
                 +-----------+-----------+
                             |
                             v
                 +-----------------------+
                 | NMR Shielding Calc    |
                 | (B3LYP/6-311+G(2d,p)) |
                 | (GIAO method)         |
                 +-----------+-----------+
                             |
                             v
                 +-----------------------+
                 |   Linear Scaling      |
                 |  (DELTA50 factors)    |
                 +-----------+-----------+
                             |
                             v
                 +-----------------------+
                 | Boltzmann Averaging   |
                 |  (Energy-weighted)    |
                 +-----------+-----------+
                             |
                             v
                 +-----------------------+
                 | Predicted Shifts      |
                 |    (1H, 13C ppm)      |
                 +-----------------------+
```

### DP4+ Context

This application implements the **chemical shift prediction** component of the DP4+ methodology. The full DP4+ analysis (developed by Goodman and Sarotti groups) uses these predictions with Bayesian statistics to assign stereochemistry. While this project does not implement the probability calculation itself, it provides the high-quality shift predictions that feed into such analyses.

---

## NMR Chemical Shift Fundamentals

### Nuclear Magnetic Resonance Basics

Certain atomic nuclei possess **nuclear spin**, making them behave like tiny magnets. When placed in an external magnetic field $B_0$, these nuclei align with or against the field and can absorb radiofrequency radiation at a characteristic **Larmor frequency**:

$$\nu_0 = \frac{\gamma B_0}{2\pi}$$

where $\gamma$ is the gyromagnetic ratio (a constant for each isotope). For $^1$H and $^{13}$C---the nuclei most commonly studied---this frequency falls in the MHz range at typical magnetic field strengths.

### Shielding and Chemical Shift

The key insight of NMR spectroscopy is that each nucleus does not experience exactly $B_0$. Instead, the electrons surrounding the nucleus create their own small magnetic field that **shields** the nucleus from the external field. The effective field felt by the nucleus is:

$$B_{eff} = B_0(1 - \sigma)$$

where $\sigma$ is the **shielding constant** (dimensionless, typically $10^{-6}$ to $10^{-4}$). Nuclei in different chemical environments have different electron densities, hence different shielding, hence different resonance frequencies.

Rather than reporting absolute frequencies (which depend on the instrument's magnetic field), we use the **chemical shift** $\delta$---the frequency difference from a reference compound, expressed in parts per million:

$$\delta = \frac{\nu_{sample} - \nu_{ref}}{\nu_{ref}} \times 10^6 \text{ ppm}$$

Since $\nu \propto (1-\sigma)$, chemical shift is related to shielding by:

$$\delta \approx \sigma_{ref} - \sigma_{sample}$$

Higher shielding (more electron density) leads to *lower* chemical shift (upfield), while lower shielding (less electron density) leads to *higher* chemical shift (downfield).

### Reference Compounds

The standard reference for $^1$H and $^{13}$C NMR is **tetramethylsilane (TMS)**, assigned $\delta = 0$ ppm for both nuclei. TMS was chosen because:

- It gives a single sharp peak (all 12 hydrogens and 4 carbons are equivalent)
- It appears upfield of most organic signals
- It is chemically inert and volatile (easily removed from samples)

In computational work, rather than calculating TMS and subtracting, we use **empirical linear scaling** (see [Linear Scaling Methodology](#linear-scaling-methodology)) which implicitly handles the reference.

### Why Quantum Mechanical Calculations?

The shielding constant depends on the complete electronic structure around each nucleus. Electrons in bonds, lone pairs, and even distant aromatic rings all contribute through:

1. **Diamagnetic shielding:** Circulation of electrons in the ground state
2. **Paramagnetic shielding:** Mixing of excited states (often deshielding)
3. **Ring current effects:** Aromatic rings create anisotropic shielding cones

Predicting these contributions requires knowing the electron density distribution, which in turn requires solving the electronic Schrodinger equation. **Density Functional Theory (DFT)** provides the most practical balance of accuracy and computational cost for molecules of interest.

**Key reference:**
Wolinski, K.; Hinton, J. F.; Pulay, P. "Efficient Implementation of the Gauge-Independent Atomic Orbital Method for NMR Chemical Shift Calculations." *J. Am. Chem. Soc.* **1990**, 112, 8251-8260.
DOI: [10.1021/ja00179a005](https://doi.org/10.1021/ja00179a005)

---

## DFT Theory Basis

### Density Functional Theory Overview

Density Functional Theory (DFT) reformulates quantum mechanics in terms of the **electron density** $\rho(\mathbf{r})$ rather than the many-electron wavefunction. The Hohenberg-Kohn theorems establish that the ground-state energy is uniquely determined by the electron density:

$$E[\rho] = T[\rho] + V_{ne}[\rho] + V_{ee}[\rho]$$

where $T$ is kinetic energy, $V_{ne}$ is nuclear-electron attraction, and $V_{ee}$ is electron-electron interaction.

The practical breakthrough came from Kohn and Sham, who introduced auxiliary orbitals to compute the kinetic energy exactly, leaving only the **exchange-correlation functional** $E_{xc}[\rho]$ to be approximated. Different choices of $E_{xc}$ define different DFT methods.

### B3LYP Hybrid Functional

This project uses **B3LYP** (Becke 3-parameter Lee-Yang-Parr), the most widely used hybrid functional for NMR calculations. It combines:

- **Exact Hartree-Fock exchange** (20%)
- **Becke gradient-corrected exchange** (72%)
- **Slater local exchange** (8%)
- **Lee-Yang-Parr correlation** (81%)
- **VWN local correlation** (19%)

The mixing parameters were fit to reproduce thermochemical data. For NMR shielding, B3LYP consistently provides accurate predictions at moderate computational cost.

**Implementation:** See [`input_gen.py`](../src/qm_nmr_calc/nwchem/input_gen.py) for how B3LYP is specified in NWChem input.

**Key reference:**
Becke, A. D. "Density-functional thermochemistry. III. The role of exact exchange." *J. Chem. Phys.* **1993**, 98, 5648-5652.
DOI: [10.1063/1.464913](https://doi.org/10.1063/1.464913)

### Basis Sets

A **basis set** is the mathematical functions used to represent molecular orbitals. Larger basis sets give more accurate results but cost more computationally. We use different basis sets for different stages:

#### Geometry Optimization: 6-31G*

For finding the minimum-energy structure, **6-31G*** provides a good balance:
- **6-31G**: Split-valence double-zeta (two functions per valence orbital)
- **\***: d-polarization functions on heavy atoms

This captures bond lengths and angles accurately enough for subsequent NMR calculations, at relatively low cost.

#### NMR Shielding: 6-311+G(2d,p)

For NMR shielding tensors, we need higher accuracy near the nuclei. **6-311+G(2d,p)** provides:
- **6-311**: Triple-zeta (three functions per valence orbital)
- **+**: Diffuse functions on heavy atoms (important for anions, lone pairs)
- **(2d,p)**: Two sets of d-polarization on heavy atoms, p-polarization on hydrogens

The diffuse functions are particularly important for magnetic properties because they better describe the electron density in the region between atoms.

**Basis set notation summary:**

| Notation | Meaning |
|----------|---------|
| 6-311 | Triple-zeta valence |
| + | Diffuse functions (sp on heavy atoms) |
| (2d,p) | Polarization: 2 d-sets on C,N,O,...; 1 p-set on H |

### GIAO Method

Computing magnetic shielding requires evaluating matrix elements involving the vector potential of the external magnetic field. A fundamental problem arises: the calculated shielding depends on where we place the **gauge origin** (the point where the vector potential is zero), but physical observables must be gauge-independent.

The **Gauge-Including Atomic Orbital (GIAO)** method solves this by attaching a field-dependent phase factor to each basis function:

$$\chi_\mu^{GIAO}(\mathbf{r}) = e^{-i\mathbf{A}_\mu \cdot \mathbf{r}} \chi_\mu(\mathbf{r})$$

where $\mathbf{A}_\mu$ is the vector potential evaluated at the center of basis function $\mu$. This makes each basis function "carry" its own gauge origin, ensuring gauge-origin independence even with finite (incomplete) basis sets.

The shielding tensor $\boldsymbol{\sigma}$ is a 3x3 matrix relating the induced field to the applied field. For NMR, we report the **isotropic shielding**:

$$\sigma_{iso} = \frac{1}{3}\text{Tr}(\boldsymbol{\sigma}) = \frac{1}{3}(\sigma_{xx} + \sigma_{yy} + \sigma_{zz})$$

This is what we convert to chemical shift via linear scaling.

**Implementation:** NWChem computes GIAO shielding via:
```
property
  shielding
end
task dft property
```

**Key reference:**
Wolinski, K.; Hinton, J. F.; Pulay, P. "Efficient Implementation of the Gauge-Independent Atomic Orbital Method for NMR Chemical Shift Calculations." *J. Am. Chem. Soc.* **1990**, 112, 8251-8260.
DOI: [10.1021/ja00179a005](https://doi.org/10.1021/ja00179a005)

---

## COSMO Solvation Model

### Implicit Solvation

NMR experiments are performed in solution, where solvent molecules affect the electronic structure of the solute. We could include explicit solvent molecules, but this would require:
- Generating representative solvent configurations
- Averaging over many configurations (expensive)
- Much larger calculations per configuration

**Implicit solvation** models instead treat the solvent as a polarizable continuum characterized by its **dielectric constant** $\varepsilon$. The solute sits in a cavity carved from this continuum. This captures the bulk electrostatic effect of solvation at minimal computational cost.

The trade-off: implicit models cannot capture specific solute-solvent interactions like hydrogen bonding. For most NMR predictions, the electrostatic effect dominates, making implicit solvation a practical choice.

### COSMO Model

The **COnductor-like Screening MOdel (COSMO)** approximates the solvent as a conductor ($\varepsilon \to \infty$) initially, then scales back to the actual dielectric:

1. Define a molecular-shaped cavity around the solute (typically van der Waals surface)
2. Place point charges on the cavity surface to screen the solute's electric field
3. In a conductor, these charges would completely screen the field
4. Scale the screening charges by factor $\frac{\varepsilon-1}{\varepsilon+0.5}$ for finite $\varepsilon$

The scaling factor approaches 1 for high $\varepsilon$ (good conductors/polar solvents) and 0 for $\varepsilon = 1$ (vacuum).

**Implementation:** See [`input_gen.py`](../src/qm_nmr_calc/nwchem/input_gen.py) for COSMO configuration:

```python
cosmo_block = f"""
cosmo
  do_gasphase False
  solvent {solvent_name}
end
"""
```

**Key reference:**
Klamt, A.; Schuurmann, G. "COSMO: A New Approach to Dielectric Screening in Solvents with Explicit Expressions for the Screening Energy and its Gradient." *J. Chem. Soc., Perkin Trans. 2* **1993**, 799-805.
DOI: [10.1039/P29930000799](https://doi.org/10.1039/P29930000799)

### Dielectric Constants

The dielectric constant determines the strength of solvation effects:

| Solvent | Dielectric Constant ($\varepsilon$) | Character |
|---------|-------------------------------------|-----------|
| Vacuum | 1.0 | No solvation |
| CHCl3 (chloroform) | 4.81 | Weakly polar |
| DMSO | 46.7 | Highly polar |
| Water | 78.4 | Highly polar |

Higher $\varepsilon$ means stronger screening of electrostatic interactions. Polar solvents stabilize polar/charged groups more, which can shift electron density and thus change NMR shielding.

**Important:** The solvent used in calculations should match experimental conditions. Our scaling factors are derived separately for each solvent (CHCl3, DMSO, vacuum) to account for solvent-specific systematic errors.

### Solvents Supported in This Project

The QM NMR Calculator supports the following solvents via NWChem's COSMO implementation:

- **chcl3** - Chloroform (CDCl3 is the common NMR solvent)
- **dmso** - Dimethyl sulfoxide (DMSO-d6)
- **water** - For aqueous NMR
- **acetone** - Acetone-d6
- **methanol** - Methanol-d4
- **vacuum** - Gas-phase calculation (no solvation)

Scaling factors are currently available for **CHCl3**, **DMSO**, and **vacuum**.

---

## Linear Scaling Methodology

### The Reference Compound Problem

The traditional approach to converting calculated shielding constants to chemical shifts uses a single reference compound (typically TMS):

$$\delta = \sigma_{ref} - \sigma_{sample}$$

However, this approach has significant limitations:

1. **TMS shielding varies:** The calculated shielding of TMS depends on the DFT method, basis set, and solvation model. A value calculated with one approach cannot be used with another.

2. **Systematic errors accumulate:** DFT methods have systematic biases that affect all calculated shieldings. A single-point reference cannot correct for these trends.

3. **Method-specific calibration needed:** Each combination of functional/basis set/solvent produces different systematic errors, requiring separate calibration.

### Empirical Regression Approach

Rather than relying on a single reference compound, **empirical linear scaling** fits a regression across many compounds with known experimental shifts. This approach:

- Uses a benchmark dataset of diverse organic molecules
- Fits method-specific scaling factors (slope and intercept)
- Implicitly handles the reference and corrects systematic errors
- Provides error statistics (MAE, RMSD) for quality assessment

The key insight is that DFT shielding errors are largely **systematic**, not random. A linear correction captures most of this systematic error.

### Mathematical Derivation

Given:
- $\sigma_{calc,i}$ = calculated shielding for atom $i$
- $\delta_{exp,i}$ = experimental chemical shift for atom $i$

We assume a linear relationship:

$$\delta_{calc} = m \cdot \sigma_{calc} + b$$

where $m$ is the slope and $b$ is the intercept. We find optimal values by minimizing the sum of squared errors (ordinary least squares):

$$\chi^2 = \sum_i \left(\delta_{exp,i} - (m \cdot \sigma_{calc,i} + b)\right)^2$$

Taking partial derivatives and setting to zero gives the **normal equations**. The solution is:

$$m = \frac{n\sum(\sigma \cdot \delta) - \sum\sigma \sum\delta}{n\sum\sigma^2 - (\sum\sigma)^2}$$

$$b = \frac{\sum\delta - m\sum\sigma}{n}$$

where $n$ is the number of data points and summations run over all benchmark compounds.

### Physical Interpretation

The scaling factors have clear physical meaning:

- **Slope near -1.0:** Indicates the DFT method is reasonably accurate. A slope of exactly -1.0 would mean no systematic scaling error.

- **Deviations from -1.0:** Reflect systematic over/underestimation by the DFT method. For example, a slope of -0.937 (as in our 1H CHCl3 factors) indicates the method slightly underestimates shielding changes.

- **Intercept:** Represents the offset from the TMS reference point. This value corresponds approximately to the calculated TMS shielding, corrected for systematic effects.

- **Example from this project:** The 1H slope of -0.9375 indicates a ~6% systematic correction. Calculated shielding differences are multiplied by 0.9375 when converting to shifts, compressing the predicted shift range slightly.

### Implementation

The linear scaling formula is implemented in [`shifts.py`](../src/qm_nmr_calc/shifts.py):

```python
# From shifts.py shielding_to_shift function
shift = slope * shielding + intercept
```

The function `get_scaling_factor()` retrieves the appropriate factors based on the DFT method, basis set, nucleus type, and solvent used in the calculation.

---

## DELTA50 Benchmark

### Dataset Description

The **DELTA50 benchmark** is a curated dataset of 50 organic molecules with diverse functional groups, used to derive empirical scaling factors for NMR chemical shift prediction.

**Key characteristics:**
- **50 molecules** spanning alcohols, ethers, ketones, amines, aromatics, and heterocycles
- **Known experimental shifts** in CDCl3 and DMSO-d6 solvents
- **~335 1H data points** and **~219 13C data points** per solvent
- **Source:** Grimblat et al. *Molecules* **2023**, 28, 2449
- **DOI:** [10.3390/molecules28062449](https://doi.org/10.3390/molecules28062449)

The dataset was designed to be representative of typical organic chemistry, avoiding problematic cases (highly strained systems, paramagnetic compounds, dynamic exchange) that would require specialized treatment.

### Derivation of Scaling Factors

The scaling factors used in this project were derived as follows:

1. **DFT calculations:** All 50 molecules computed with B3LYP/6-311+G(2d,p) and COSMO solvation (CHCl3 or DMSO)

2. **Linear regression:** Calculated shieldings regressed against experimental shifts for each nucleus type (1H, 13C) and each solvent

3. **Outlier removal:** 3-sigma criterion applied to remove significant outliers (7 for 1H, 2 for 13C in CHCl3). These typically represent problematic functional groups or experimental errors.

4. **Quality metrics:** R-squared, MAE, and RMSD computed on the remaining data

### Project Scaling Factors

The following scaling factors are used in this project, stored in [`scaling_factors.json`](../src/qm_nmr_calc/data/scaling_factors.json):

| Parameter | 1H (CHCl3) | 13C (CHCl3) | 1H (DMSO) | 13C (DMSO) | 1H (vacuum) | 13C (vacuum) |
|-----------|------------|-------------|-----------|------------|-------------|--------------|
| Slope | -0.9375 | -0.9497 | -0.9323 | -0.9429 | -0.9554 | -0.9726 |
| Intercept | 29.92 | 172.69 | 29.73 | 171.77 | 30.54 | 175.71 |
| R^2 | 0.9952 | 0.9978 | 0.9951 | 0.9974 | 0.9934 | 0.9980 |
| MAE | 0.12 ppm | 1.95 ppm | 0.13 ppm | 2.15 ppm | 0.15 ppm | 1.74 ppm |
| RMSD | 0.16 ppm | 2.69 ppm | 0.17 ppm | 2.92 ppm | 0.19 ppm | 2.56 ppm |
| N points | 335 | 219 | 335 | 219 | 336 | 219 |
| Outliers removed | 7 | 2 | 7 | 2 | 6 | 2 |

**Interpreting the table:**

- **High R-squared (>0.99):** Indicates excellent linear correlation between calculated and experimental values
- **Low MAE:** Mean absolute error represents typical prediction accuracy (0.12-0.15 ppm for 1H, 1.7-2.2 ppm for 13C)
- **RMSD > MAE:** Indicates some larger errors exist, but the distribution is not highly skewed

### Comparison Across Solvents

The scaling factors vary slightly between solvents:

- **Vacuum vs solvated:** Vacuum calculations show steeper slopes (closer to -1.0) and larger intercepts, reflecting the absence of solvation screening
- **CHCl3 vs DMSO:** Small differences (<1%) in slope, with DMSO showing slightly larger MAE (expected for the more polar solvent with stronger solute-solvent interactions)

**Important:** Always use scaling factors that match your calculation conditions. Using CHCl3 factors for a DMSO calculation will introduce systematic errors.
