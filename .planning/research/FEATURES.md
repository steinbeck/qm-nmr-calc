# Feature Landscape: NMReData Export for Predicted NMR Data

**Domain:** NMR chemical shift prediction with computational chemistry
**Researched:** 2026-02-01
**Confidence:** HIGH (official NMReData documentation, format v1.0/1.1)

## Executive Summary

NMReData is a standard format for associating NMR parameters with molecular structures using SDF (Structure Data Format) with custom tags. The format explicitly supports **both experimental and computational/predicted data**, making it well-suited for qm-nmr-calc's DFT-predicted chemical shifts.

For predicted data, NMReData requires:
1. **Core structure data** (MOL block with 3D coordinates)
2. **Metadata tags** (version, solvent, temperature)
3. **Assignment mapping** (atom indices to chemical shifts)
4. **Optional spectrum-like representation** (even though no actual spectrum exists)

The format is flexible - many experimental-only fields (coupling constants, peak integrals, multiplicities) are optional and can be omitted for purely computational data.

## Table Stakes Fields

Required fields for a **minimal valid NMReData file** with predicted data.

| Field | Purpose | Data Mapping | Format | Notes |
|-------|---------|--------------|--------|-------|
| **MOL block** | 3D molecular structure | From optimized geometry XYZ | V2000/V3000 SDF | WE HAVE: `optimized.xyz` from DFT |
| `NMREDATA_VERSION` | Format version identifier | Static: "1.1" | Single line: `1.1\` | Current stable version |
| `NMREDATA_LEVEL` | Assignment ambiguity level | Static: "0" (no ambiguity) | Single line: `0` | We have full atom assignments |
| `NMREDATA_SOLVENT` | Solvent for calculation | From job params: `chcl3` → `CDCl3` | e.g., `CDCl3\` | Map our codes to NMR solvents |
| `NMREDATA_TEMPERATURE` | Calculation temperature | Static or from ensemble: `298.15 K` | e.g., `298.15 K\` | Default 298.15 K for DFT |
| `NMREDATA_ASSIGNMENT` | Atom index → shift mapping | From `h1_shifts`, `c13_shifts` | `label, shift, atom_idx\` | **Core data export** |
| `NMREDATA_FORMULA` | Molecular formula | From RDKit molecule | e.g., `C2H6O` | Optional but strongly recommended |
| `NMREDATA_SMILES` | SMILES representation | From input SMILES | Isomeric preferred | Optional but strongly recommended |

**Complexity:** LOW - straightforward mapping from existing data structures

## Optional But Valuable Fields

Fields that enhance usefulness and interoperability but aren't strictly required.

| Field | Purpose | Data Mapping | Value | Why Include |
|-------|---------|--------------|-------|-------------|
| `NMREDATA_ID` | Metadata and provenance | Custom fields | See below | Track calculation origin |
| `NMREDATA_INCHI` | Structure identifier | From RDKit | InChI string | Better than SMILES for matching |
| `NMREDATA_1D_1H` | 1H spectrum representation | From `h1_shifts` | Pseudo-spectrum | Enable visualization tools |
| `NMREDATA_1D_13C` | 13C spectrum representation | From `c13_shifts` | Pseudo-spectrum | Enable visualization tools |
| `NMREDATA_J` | Coupling constants | NOT AVAILABLE | Omit | We don't calculate J-couplings |

### NMREDATA_ID Subfields (Recommended)

The `NMREDATA_ID` tag supports custom key-value pairs for provenance tracking:

| Subfield | Value | Purpose |
|----------|-------|---------|
| `Title` | From `input_name` or auto-generated | Human-readable identifier |
| `Source` | "Computational" | Distinguish from experimental |
| `Software` | "qm-nmr-calc vX.X.X" | Track generating software |
| `Method` | "DFT/GIAO/B3LYP" | Calculation method |
| `Basis_Geometry` | "6-31G*" or "6-311+G(2d,p)" | From preset |
| `Basis_NMR` | "6-31G*" or "6-311+G(2d,p)" | From preset |
| `Scaling` | "DELTA50-derived OLS" | Scaling method used |
| `Job_ID` | Job ID from qm-nmr-calc | Traceability back to calculation |
| `Comment` | "Predicted chemical shifts" | Clarify data nature |

**Example:**
```
>  <NMREDATA_ID>
Title=Ethanol predicted NMR\
Source=Computational\
Software=qm-nmr-calc v2.3.0\
Method=DFT/GIAO/B3LYP/COSMO\
Basis_Geometry=6-311+G(2d,p)\
Basis_NMR=6-311+G(2d,p)\
Scaling=DELTA50-derived OLS\
Job_ID=a1b2c3d4e5f6\
Comment=Predicted 1H/13C chemical shifts from quantum chemistry\
```

**Complexity:** LOW - string formatting from metadata

### Pseudo-Spectrum Fields

Even though we don't have actual spectra, we can generate **pseudo-spectrum tags** to enable visualization in NMReData-compatible tools:

**NMREDATA_1D_1H:**
```
>  <NMREDATA_1D_1H>
Larmor=400.000000\
Spectrum_Location=none\
7.2345, L=H1\
4.5678, L=H2\
1.2345, L=H3\
```

**NMREDATA_1D_13C:**
```
>  <NMREDATA_1D_13C>
Larmor=100.000000\
Spectrum_Location=none\
128.45, L=C1\
65.32, L=C2\
18.77, L=C3\
```

- `Larmor`: Arbitrary but realistic frequency (400 MHz for 1H, 100 MHz for 13C)
- `Spectrum_Location=none`: No actual spectrum file
- Signal lines: `shift, L=label` (minimal format, no multiplicity/integrals)

**Value:** Enables tools like Mestrelab Mnova to display predicted shifts visually

**Complexity:** MEDIUM - requires generating labels and formatting per-atom data

## Not Applicable Fields

Experimental-only fields we **skip entirely** for predicted data.

| Field | Why Skip | Alternative |
|-------|----------|-------------|
| `NMREDATA_J` | No J-coupling predictions | Omit - not part of our calculation |
| `NMREDATA_2D_*` | No 2D spectra | Omit - experimental technique only |
| Multiplicity (`S=`) | No peak shapes from DFT | Omit from 1D tags |
| Integral (`E=`, `N=`) | No experimental intensities | Omit from 1D tags |
| Peak width (`W=`) | No lineshape data | Omit from 1D tags |
| `T1`, `T2`, `Diff` | No relaxation/diffusion data | Omit - experimental measurements |
| `Spectrum_Location` (with file) | No FID/spectra files | Use `Spectrum_Location=none` |
| `MD5_*` | No spectrum files to hash | Omit checksum fields |
| `Pulseprogram` | No pulse sequence | Omit from 1D tags |
| `CorType` | No 2D experiments | N/A |

**Note:** The NMReData format is designed to be flexible. Omitting experimental-only fields is **explicitly supported** for computational data sources.

## Data We Have → NMReData Field Mapping

Comprehensive mapping of qm-nmr-calc outputs to NMReData fields.

### From Job Metadata (`metadata.json`)

| Our Data | NMReData Field | Transformation |
|----------|---------------|----------------|
| `input_smiles` | `NMREDATA_SMILES` | Direct copy (isomeric SMILES preferred) |
| `input_name` | `NMREDATA_ID` → `Title` | Use if provided, else generate from formula |
| `solvent` | `NMREDATA_SOLVENT` | Map: `chcl3`→`CDCl3`, `dmso`→`DMSO-d6`, `vacuum`→`none` |
| `preset` | `NMREDATA_ID` → `Basis_*` | Map: `draft`→`6-31G*`, `production`→`6-311+G(2d,p)` |
| `functional` | `NMREDATA_ID` → `Method` | Static: `B3LYP` |
| `job_id` | `NMREDATA_ID` → `Job_ID` | Direct copy |

### From NMR Results (`results.json` / API)

| Our Data | NMReData Field | Transformation |
|----------|---------------|----------------|
| `h1_shifts[].index` | `NMREDATA_ASSIGNMENT` | Atom number (1-based, already correct) |
| `h1_shifts[].shift` | `NMREDATA_ASSIGNMENT` | Chemical shift in ppm (4+ decimals) |
| `h1_shifts[].atom` | Label generation | Generate label: `H{index}` |
| `c13_shifts[].index` | `NMREDATA_ASSIGNMENT` | Atom number (1-based) |
| `c13_shifts[].shift` | `NMREDATA_ASSIGNMENT` | Chemical shift in ppm (4+ decimals) |
| `c13_shifts[].atom` | Label generation | Generate label: `C{index}` |

### From Geometry Files

| Our Data | NMReData Field | Transformation |
|----------|---------------|----------------|
| `optimized.xyz` | MOL block in SDF | Convert XYZ → MOL format (RDKit can do this) |
| Atom coordinates | MOL block V2000 | 3D coordinates from optimized geometry |

### From RDKit Molecule

| Our Data | NMReData Field | Method |
|----------|---------------|--------|
| Molecular formula | `NMREDATA_FORMULA` | `rdkit.Chem.rdMolDescriptors.CalcMolFormula()` |
| InChI | `NMREDATA_INCHI` | `rdkit.Chem.inchi.MolToInchi()` |
| SMILES (canonical) | `NMREDATA_SMILES` | `rdkit.Chem.MolToSmiles()` |

### From Ensemble Metadata (if ensemble mode)

| Our Data | NMReData Field | Use |
|----------|---------------|-----|
| `temperature_k` | `NMREDATA_TEMPERATURE` | Use ensemble temperature instead of 298.15 K |
| `conformer_count` | `NMREDATA_ID` → `Comment` | Mention in comment: "Boltzmann-weighted average of N conformers" |
| `method` (rdkit_kdg/crest) | `NMREDATA_ID` → `Comment` | Mention conformer generation method |

## Ensemble vs Single-Conformer Handling

| Scenario | Temperature | Comment Field |
|----------|-------------|---------------|
| **Single conformer** | `298.15 K` (default) | "Predicted from single-conformer calculation" |
| **Ensemble (RDKit)** | From `ensemble_metadata.temperature_k` | "Boltzmann-weighted average of {N} RDKit KDG conformers" |
| **Ensemble (CREST)** | From `ensemble_metadata.temperature_k` | "Boltzmann-weighted average of {N} CREST/GFN2-xTB conformers" |

## Feature Dependencies

```
Required dependencies (must implement in order):
1. XYZ → MOL conversion ────┐
2. SMILES/formula extraction│
3. ASSIGNMENT tag generator │
                            ├──> Minimal valid NMReData file
4. VERSION/LEVEL/SOLVENT/   │
   TEMPERATURE generators   │
5. SDF file assembly ────────┘

Optional enhancements (add value, not required):
6. NMREDATA_ID with provenance → Better traceability
7. NMREDATA_INCHI → Better structure matching
8. 1D pseudo-spectrum tags → Visualization tool compatibility
```

## MVP Recommendation

**For MVP (minimal valid NMReData export):**

1. **MOL block** from `optimized.xyz` (RDKit: XYZ → MOL)
2. **Required tags:**
   - `NMREDATA_VERSION` = "1.1"
   - `NMREDATA_LEVEL` = "0"
   - `NMREDATA_SOLVENT` (mapped from job solvent)
   - `NMREDATA_TEMPERATURE` = "298.15 K" (or ensemble temp)
   - `NMREDATA_ASSIGNMENT` (from `h1_shifts` + `c13_shifts`)
3. **Strongly recommended:**
   - `NMREDATA_FORMULA` (from RDKit)
   - `NMREDATA_SMILES` (from input or RDKit)

**Defer to post-MVP:**
- `NMREDATA_ID` with full provenance metadata (adds value but not required for validation)
- `NMREDATA_INCHI` (nice-to-have for structure matching)
- `NMREDATA_1D_1H` / `NMREDATA_1D_13C` pseudo-spectrum tags (for tool compatibility)

**Rationale:** MVP focuses on core assignment export. Provenance and visualization enhancements can be added incrementally without breaking existing files.

## Label Generation Strategy

NMReData requires labels in the `NMREDATA_ASSIGNMENT` tag. We need a labeling scheme for atoms.

### Recommended Approach: Simple Index-Based Labels

```
H atoms: H1, H2, H3, ... (based on atom index)
C atoms: C1, C2, C3, ... (based on atom index)
```

**Example for ethanol (CH3CH2OH):**
```
>  <NMREDATA_ASSIGNMENT>
C1, 18.7704, 1\
C2, 63.5132, 2\
H1, 1.2436, 3\
H2, 1.2436, 4\
H3, 1.2436, 5\
H4, 3.8300, 6\
H5, 3.8300, 7\
O1, 101.4098, 8\
H6, 0.3412, 9\
```

**Alternative: Chemical environment labels** (e.g., `CH3`, `CH2`, `OH`)
- More readable but harder to generate algorithmically
- Requires chemical interpretation (which CH3 is which?)
- Defer to future enhancement

**MVP decision:** Use index-based labels (`H{index}`, `C{index}`). Simple, unambiguous, maps directly to atom indices.

## Solvent Code Mapping

| qm-nmr-calc Code | NMReData Standard Name | Notes |
|------------------|------------------------|-------|
| `chcl3` | `CDCl3` | Standard chloroform-d |
| `dmso` | `DMSO-d6` | Deuterated DMSO |
| `vacuum` | `none` | Gas phase, no solvent |

**Future solvents:** When added to qm-nmr-calc, map to deuterated NMR equivalents (e.g., `acetone` → `Acetone-d6`).

## Format Precision Requirements

| Field Type | Precision | Example |
|------------|-----------|---------|
| Chemical shifts | **4+ decimals** | `7.2345` (not `7.23`) |
| Temperature | 2 decimals + unit | `298.15 K` |
| Atom indices | Integer | `1, 2, 3` |
| Labels | Alphanumeric | `H1`, `C2` |

**Critical:** NMReData specification **requires 4 decimal places minimum** for chemical shifts. Our current data has this precision already.

## Validation Criteria

A valid NMReData file for qm-nmr-calc must:

1. ✅ Parse as valid SDF (MOL block + tags)
2. ✅ Contain `NMREDATA_VERSION` = "1.1"
3. ✅ Contain `NMREDATA_LEVEL` = "0"
4. ✅ Contain `NMREDATA_ASSIGNMENT` with all predicted atoms
5. ✅ Have chemical shifts with ≥4 decimal places
6. ✅ Reference atom indices matching MOL block atom count
7. ✅ Include `NMREDATA_SOLVENT` and `NMREDATA_TEMPERATURE`

**Testing strategy:** Use NMReData-compatible tools (Mestrelab Mnova, nmrshiftdb2) to verify files load correctly.

## Tool Compatibility Targets

| Tool | Purpose | NMReData Support | Priority |
|------|---------|------------------|----------|
| **Mestrelab Mnova** | NMR analysis software | Full support (v1.0/1.1) | High |
| **nmrshiftdb2** | Open NMR database | Full support | High |
| **ACD/Labs NMR** | Commercial NMR software | Full support | Medium |
| **ChemDraw** | Structure drawing | SDF import (partial) | Low |

**Goal:** Files should load in Mnova and nmrshiftdb2 without errors. Full visualization may require pseudo-spectrum tags (post-MVP).

## Sources

### Primary Documentation
- [NMReData Tag Format (v1.0/1.1)](https://github.com/NMReDATAInitiative/NMReDATAInitiative.github.io/blob/master/NMReDATA_tag_format.md) - Official format specification
- [NMReData 1D Attributes](https://github.com/NMReDATAInitiative/NMReDATAInitiative.github.io/blob/master/1D_attributes.md) - Spectrum signal attributes
- [NMReData Examples Repository](https://github.com/NMReDATAInitiative/Examples-of-NMR-records) - Example files including computational data

### Literature
- [NMReDATA, a standard to report the NMR assignment and parameters of organic compounds](https://pmc.ncbi.nlm.nih.gov/articles/PMC6226248/) - Primary publication (Magnetic Resonance in Chemistry, 2018)
- [NMReData Initiative Website](https://nmredata.org/) - Official initiative homepage (currently offline)

### Computational NMR Context
- [Accurate Prediction of NMR Chemical Shifts: Integrating DFT Calculations with 3D GNNs](https://pmc.ncbi.nlm.nih.gov/articles/PMC11209944/) - DFT-GIAO methodology
- [IR-NMR multimodal computational spectra dataset](https://www.nature.com/articles/s41597-025-05729-8) - Computational NMR dataset (2025)
- [nmrshiftdb2 - open nmr database](https://nmrshiftdb.nmr.uni-koeln.de/) - Supports predicted/calculated NMR data

**Confidence level:** HIGH - Official documentation accessed, format specification version 1.0/1.1 is stable and widely adopted. Example files confirm computational data is explicitly supported.
