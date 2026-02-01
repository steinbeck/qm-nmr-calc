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
