# Phase 8: DELTA50 Setup - Research

**Researched:** 2026-01-21
**Domain:** Computational chemistry benchmarking, NMR shift prediction validation
**Confidence:** HIGH

## Summary

Phase 8 establishes the DELTA50 benchmark infrastructure, which consists of 50 carefully curated small organic molecules with experimental ¹H and ¹³C NMR chemical shifts. The benchmark was published in March 2023 in *Molecules* (DOI: 10.3390/molecules28062449) and provides supporting information including XYZ coordinate files and experimental shifts in Excel format.

The standard approach for benchmark infrastructure combines:
1. **Data repository**: Commit DELTA50 structures and experimental shifts to the codebase (not runtime downloads)
2. **Python CLI with argparse**: Subcommand pattern (`python -m qm_nmr_calc.benchmark run`)
3. **RDKit for file I/O**: Parse XYZ/SDF files, write NMReData-format SDF outputs
4. **Pydantic models**: Type-safe data structures for molecules and results
5. **Resume-capable runner**: Check for existing outputs, skip completed calculations

**Primary recommendation:** Use RDKit's native SDF reading/writing with custom property tags for NMR data. Build a sequential benchmark runner with tqdm progress tracking, JSONL intermediate results, and final CSV/SDF outputs. Implement resume capability by checking for existing output directories before running each molecule×method×solvent combination.

## Standard Stack

The established libraries/tools for chemistry benchmarking in Python:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| RDKit | 2025.9.3+ | Molecule I/O, structure handling | De facto standard for cheminformatics in Python, already in project dependencies |
| Pydantic | 2.12.5+ | Data validation, serialization | Already used in project for JobInput/JobStatus models, type-safe JSON handling |
| argparse | stdlib | CLI argument parsing | Standard library, zero dependencies, excellent subcommand support |
| tqdm | Latest | Progress bars | Industry standard for CLI progress tracking, ~60ns overhead per iteration |
| pandas | Latest | CSV/DataFrame output | Standard for tabular scientific data in Python |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| pathlib | stdlib | Path manipulation | All file operations (mkdir, exists, glob) |
| json | stdlib | JSON serialization | Pydantic model serialization, JSONL writing |
| logging | stdlib | Structured logging | Failure tracking, debug information |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| argparse | click or typer | Click/typer offer cleaner APIs with decorators, but argparse requires zero dependencies and has excellent subcommand support for git-like CLIs |
| CSV output | JSON only | CSV is human-readable in Excel, essential for chemistry workflows |
| RDKit SDF | OpenBabel | OpenBabel can convert formats but doesn't preserve molecule objects for property manipulation; RDKit already in dependencies |

**Installation:**
```bash
# Core dependencies already in pyproject.toml
# Add for benchmark runner:
uv add tqdm pandas
```

## Architecture Patterns

### Recommended Project Structure
```
src/qm_nmr_calc/
├── benchmark/
│   ├── __init__.py         # Exports run_benchmark()
│   ├── __main__.py         # CLI entry point: python -m qm_nmr_calc.benchmark
│   ├── data_loader.py      # Load DELTA50 molecules and experimental shifts
│   ├── models.py           # Pydantic models for benchmark data
│   └── runner.py           # Sequential calculation runner with resume
data/
├── benchmark/
│   ├── delta50/
│   │   ├── molecules/      # 50 XYZ files (from supporting info)
│   │   └── experimental_shifts.json  # Parsed from Excel SI
│   └── results/
│       ├── molecule_01/
│       │   ├── B3LYP_CHCl3/
│       │   │   ├── input.nw
│       │   │   ├── output.out
│       │   │   └── shifts.json
│       │   └── WP04_DMSO/
│       └── summary.csv     # Aggregated results
```

### Pattern 1: Subcommand CLI with argparse
**What:** Git-like CLI with `run`, `status`, `clean` subcommands
**When to use:** Multi-action CLIs where each subcommand has different arguments
**Example:**
```python
# Source: Python docs - https://docs.python.org/3/library/argparse.html
import argparse

def main():
    parser = argparse.ArgumentParser(prog="qm_nmr_calc.benchmark")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # run subcommand
    run_parser = subparsers.add_parser("run", help="Execute benchmark calculations")
    run_parser.add_argument("--resume", action="store_true", help="Skip completed calculations")
    run_parser.add_argument("--molecules", nargs="+", help="Specific molecule IDs to run")

    # status subcommand
    status_parser = subparsers.add_parser("status", help="Show benchmark progress")

    args = parser.parse_args()

    if args.command == "run":
        run_benchmark(resume=args.resume, molecules=args.molecules)
    elif args.command == "status":
        show_status()

if __name__ == "__main__":
    main()
```

### Pattern 2: Resume-Capable Batch Processing
**What:** Check for existing outputs before executing, track progress in JSONL
**When to use:** Long-running calculations that may fail mid-run
**Example:**
```python
from pathlib import Path
from tqdm import tqdm

def run_benchmark(resume: bool = True):
    results_dir = Path("data/benchmark/results")
    progress_file = results_dir / "progress.jsonl"

    # Build task list: 50 molecules × 2 functionals × 2 solvents = 200 calculations
    tasks = build_task_matrix()

    if resume:
        completed = load_completed_from_jsonl(progress_file)
        tasks = [t for t in tasks if t["id"] not in completed]

    for task in tqdm(tasks, desc="Running benchmark"):
        output_dir = results_dir / task["molecule"] / f"{task['functional']}_{task['solvent']}"

        if output_dir.exists() and (output_dir / "shifts.json").exists():
            continue  # Already completed

        try:
            result = run_calculation(task)
            save_result(output_dir, result)
            append_to_jsonl(progress_file, {"id": task["id"], "status": "complete"})
        except Exception as e:
            logging.error(f"Failed {task['id']}: {e}")
            append_to_jsonl(progress_file, {"id": task["id"], "status": "failed", "error": str(e)})
```

### Pattern 3: Pydantic Models for Benchmark Data
**What:** Type-safe data structures for molecules, experimental shifts, and results
**When to use:** Always - provides validation, serialization, and documentation
**Example:**
```python
from pydantic import BaseModel

class MoleculeData(BaseModel):
    """DELTA50 molecule with experimental shifts."""
    id: str  # e.g., "molecule_01"
    name: str  # IUPAC or common name
    xyz_file: str  # Relative path to XYZ file
    h1_shifts: list[float]  # Experimental 1H shifts in ppm
    c13_shifts: list[float]  # Experimental 13C shifts in ppm

class BenchmarkResult(BaseModel):
    """Single calculation result."""
    molecule_id: str
    functional: str
    basis_set: str
    solvent: str
    calculated_h1: list[float]
    calculated_c13: list[float]
    h1_mae: float  # Mean absolute error vs experimental
    c13_mae: float

# Serialize to JSON
result.model_dump_json()

# Write to JSONL
with open("progress.jsonl", "a") as f:
    f.write(result.model_dump_json() + "\n")
```

### Pattern 4: SDF Output with Custom Properties
**What:** Write calculation results to SDF files with NMR shift properties
**When to use:** Chemistry-friendly output for visualization in ChemDraw/Maestro
**Example:**
```python
from rdkit import Chem

# Load molecule from XYZ
mol, xyz_block = load_geometry_file("molecule_01.xyz")

# Add NMR shifts as properties
for i, shift in enumerate(h1_shifts):
    mol.SetProp(f"NMR_1H_{i+1}", f"{shift:.2f}")

for i, shift in enumerate(c13_shifts):
    mol.SetProp(f"NMR_13C_{i+1}", f"{shift:.2f}")

# Add metadata
mol.SetProp("Functional", "B3LYP")
mol.SetProp("Solvent", "CHCl3")
mol.SetProp("BasisSet", "6-311+G(2d,p)")

# Write SDF with properties
with Chem.SDWriter("molecule_01_B3LYP_CHCl3.sdf") as w:
    w.write(mol)
```

### Anti-Patterns to Avoid
- **Don't re-download data at runtime**: Commit DELTA50 structures/shifts to repo for reproducibility
- **Don't use shell globbing for completion checks**: Use explicit output file existence checks (e.g., `shifts.json` presence)
- **Don't store all results in memory**: Write incrementally to JSONL, aggregate at end
- **Don't use relative paths**: Always use absolute paths or Path.cwd()-relative for portability

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Molecule file parsing | Custom XYZ parser | RDKit's `MolFromXYZFile()`, `MolFromMolFile()` | Handles edge cases (charges, radicals), integrates with rest of RDKit ecosystem |
| Progress tracking | Print statements with counters | tqdm with context manager | Automatic ETA calculation, no terminal pollution, works with logging |
| CSV writing | Manual string formatting | pandas `DataFrame.to_csv()` | Handles escaping, type conversion, append mode for incremental writes |
| Path operations | os.path string manipulation | pathlib.Path | Type-safe, cross-platform, cleaner API (`mkdir(parents=True, exist_ok=True)`) |
| JSON serialization | Manual dict construction | Pydantic `model_dump_json()` | Type validation, automatic datetime handling, schema generation |
| CLI argument validation | Manual if/else checks | argparse choices and required flags | Automatic help text, error messages, type coercion |

**Key insight:** RDKit already provides comprehensive molecule I/O. Don't parse chemical file formats manually - use RDKit's battle-tested parsers and let it handle bond perception, formal charges, and coordinate validation.

## Common Pitfalls

### Pitfall 1: XYZ Files Require Bond Perception
**What goes wrong:** Loading XYZ files without bond determination leads to molecules with no bonds, breaking property calculations
**Why it happens:** XYZ format only contains coordinates, not connectivity. RDKit needs explicit bond determination.
**How to avoid:** Always call `rdDetermineBonds.DetermineBonds(mol, charge=0)` after loading XYZ files
**Warning signs:** Molecule displays as disconnected atoms, `GetNumBonds()` returns 0

### Pitfall 2: DELTA50 Experimental Shifts Are Atom-Assignment Specific
**What goes wrong:** Comparing calculated vs experimental shifts without proper atom mapping gives meaningless errors
**Why it happens:** DELTA50 provides assigned chemical shifts, but atom ordering in XYZ files may differ from assignment order
**How to avoid:** DELTA50 supporting info includes fully assigned NMR spectra - use these assignments to map calculated atoms to experimental shifts. Store mapping explicitly in `experimental_shifts.json`.
**Warning signs:** Unreasonably high MAE values (>10 ppm for 1H, >50 ppm for 13C)

### Pitfall 3: Resume Logic Must Check All Output Files
**What goes wrong:** Checking only for output directory existence causes skipping failed calculations that created dirs but no results
**Why it happens:** NWChem creates scratch directories even if calculation fails
**How to avoid:** Check for specific marker file (`shifts.json` or final SDF file) in addition to directory existence
**Warning signs:** "Resume" shows 100% complete but no actual results files exist

### Pitfall 4: NMReData Format Is Not Automatic
**What goes wrong:** Expecting SDF files to automatically contain NMR data in standard format
**Why it happens:** NMReData is a convention (tags prefixed with `NMREDATA_`), not automatically handled by RDKit
**How to avoid:** Manually set molecule properties with `NMREDATA_` prefix if targeting NMReData compatibility, or use simple custom tags for internal use
**Warning signs:** Chemistry software doesn't recognize NMR data in SDF files

### Pitfall 5: Benchmark Molecule IDs Must Match NWChem Job Names
**What goes wrong:** Using molecule names with spaces/special chars as NWChem job names causes parsing failures
**Why it happens:** NWChem expects job names to be valid identifiers (no spaces, special chars)
**How to avoid:** Use sanitized IDs (`molecule_01`) for NWChem jobs, store human-readable names separately in metadata
**Warning signs:** NWChem outputs parse correctly but can't match to input molecules

## Code Examples

Verified patterns from official sources and current codebase:

### Loading DELTA50 Molecule Structures
```python
# Source: Current codebase - src/qm_nmr_calc/nwchem/geometry.py
from pathlib import Path
from qm_nmr_calc.nwchem import load_geometry_file

def load_delta50_molecules(data_dir: Path) -> dict[str, tuple]:
    """Load all 50 DELTA50 molecules from XYZ files.

    Returns:
        Dict mapping molecule_id to (RDKit Mol, XYZ block) tuple
    """
    molecules = {}
    xyz_dir = data_dir / "benchmark" / "delta50" / "molecules"

    for xyz_file in sorted(xyz_dir.glob("*.xyz")):
        molecule_id = xyz_file.stem  # e.g., "molecule_01"
        mol, xyz_block = load_geometry_file(xyz_file, charge=0)
        molecules[molecule_id] = (mol, xyz_block)

    return molecules
```

### Incremental JSONL Result Writing
```python
# Source: JSON Lines spec - https://jsonlines.org/
from pathlib import Path
import json

def append_to_jsonl(filepath: Path, data: dict):
    """Append a result to JSONL progress file."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with filepath.open("a") as f:
        f.write(json.dumps(data) + "\n")

def load_completed_from_jsonl(filepath: Path) -> set[str]:
    """Load completed task IDs from JSONL file."""
    if not filepath.exists():
        return set()

    completed = set()
    with filepath.open("r") as f:
        for line in f:
            data = json.loads(line)
            if data.get("status") == "complete":
                completed.add(data["id"])

    return completed
```

### Calculating Mean Absolute Error
```python
# Standard pattern for benchmark validation
import numpy as np

def calculate_mae(calculated: list[float], experimental: list[float]) -> float:
    """Calculate mean absolute error between calculated and experimental shifts.

    Args:
        calculated: Calculated chemical shifts in ppm
        experimental: Experimental chemical shifts in ppm (must be same order)

    Returns:
        Mean absolute error in ppm

    Raises:
        ValueError: If lists have different lengths
    """
    if len(calculated) != len(experimental):
        raise ValueError(
            f"Length mismatch: {len(calculated)} calculated vs {len(experimental)} experimental"
        )

    return float(np.mean(np.abs(np.array(calculated) - np.array(experimental))))
```

### Directory Structure Creation
```python
# Source: pathlib documentation - https://docs.python.org/3/library/pathlib.html
from pathlib import Path

def setup_benchmark_directories(base_dir: Path, molecule_ids: list[str],
                                functionals: list[str], solvents: list[str]):
    """Create directory structure for benchmark results.

    Creates: data/benchmark/results/molecule_XX/FUNCTIONAL_SOLVENT/
    """
    results_dir = base_dir / "benchmark" / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    for mol_id in molecule_ids:
        mol_dir = results_dir / mol_id
        mol_dir.mkdir(exist_ok=True)

        for functional in functionals:
            for solvent in solvents:
                method_dir = mol_dir / f"{functional}_{solvent}"
                method_dir.mkdir(exist_ok=True)
```

### Aggregating Results to CSV
```python
# Source: pandas documentation - https://pandas.pydata.org/docs/
import pandas as pd
from pathlib import Path

def aggregate_results_to_csv(results_dir: Path, output_file: Path):
    """Aggregate all benchmark results into summary CSV."""
    records = []

    for mol_dir in sorted(results_dir.glob("molecule_*")):
        for method_dir in mol_dir.glob("*_*"):
            shifts_file = method_dir / "shifts.json"
            if not shifts_file.exists():
                continue

            with shifts_file.open() as f:
                data = json.load(f)

            # Parse directory names
            molecule_id = mol_dir.name
            functional, solvent = method_dir.name.split("_")

            records.append({
                "molecule_id": molecule_id,
                "functional": functional,
                "solvent": solvent,
                "h1_mae": data.get("h1_mae"),
                "c13_mae": data.get("c13_mae"),
                "num_h1": len(data.get("h1_shifts", [])),
                "num_c13": len(data.get("c13_shifts", [])),
            })

    df = pd.DataFrame(records)
    df.to_csv(output_file, index=False)
    print(f"Summary written to {output_file} ({len(records)} results)")
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Download benchmark data at runtime | Commit to repository | Ongoing trend (2023+) | Reproducibility: ensures exact dataset version used |
| Single JSON blob for results | JSONL streaming writes | ~2020-2025 | Memory efficiency, better failure recovery |
| os.path string manipulation | pathlib.Path objects | Python 3.4+ (2014) | Type safety, cross-platform compatibility |
| Manual CLI parsing | argparse with subparsers | Python 3.7+ subparsers.required | Cleaner git-like CLI patterns |
| TMS absolute referencing | Empirical scaling factors (Pierens 2014) | 2014+ | Accuracy: accounts for systematic DFT errors |

**Deprecated/outdated:**
- **OpenBabel Python bindings**: RDKit has become the standard for Python cheminformatics (2020+)
- **Pickle for model serialization**: Pydantic JSON serialization is safer and more portable (2021+)
- **argparse.FileType**: Use pathlib.Path instead for better error handling and type safety (2022+)

## Open Questions

Things that couldn't be fully resolved:

1. **DELTA50 Atom Assignment Mapping**
   - What we know: Supporting info includes "fully assigned experimental NMR spectra"
   - What's unclear: Exact format of assignment data (PDF spectra vs machine-readable mapping)
   - Recommendation: Download SI Excel file, examine structure. If assignments aren't atom-indexed, use chemical shift proximity matching (nearest calculated to experimental within 0.5 ppm for 1H, 5 ppm for 13C)

2. **WP04 Functional Scaling Factors**
   - What we know: DELTA50 paper found WP04/6-311++G(2d,p) best for 1H shifts
   - What's unclear: User wants both B3LYP and WP04 benchmarked, but current codebase only has B3LYP scaling factors
   - Recommendation: Check DELTA50 SI Excel file for WP04 scaling factors. If not present, use raw TMS referencing (sigma_TMS - sigma_calc) for WP04 initially, note in results that empirical scaling could improve accuracy.

3. **Experimental Solvent Matching**
   - What we know: DELTA50 measured in CDCl₃ at 298K
   - What's unclear: Should benchmark run DMSO calculations when experimental is CDCl₃-only?
   - Recommendation: Run both CHCl₃ and DMSO as planned (tests method transferability), but MAE calculations only valid for CHCl₃ vs experimental. Document DMSO as "calculated only, no experimental reference" in outputs.

## Sources

### Primary (HIGH confidence)
- DELTA50 publication: [https://pmc.ncbi.nlm.nih.gov/articles/PMC10051451/](https://pmc.ncbi.nlm.nih.gov/articles/PMC10051451/) - Dataset composition, experimental conditions
- DELTA50 supporting info: [https://www.mdpi.com/article/10.3390/molecules28062449/s1](https://www.mdpi.com/article/10.3390/molecules28062449/s1) - XYZ files, experimental shifts
- RDKit 2025.09.4 documentation: [https://www.rdkit.org/docs/GettingStartedInPython.html](https://www.rdkit.org/docs/GettingStartedInPython.html) - SDF parsing, molecule I/O
- Python pathlib docs: [https://docs.python.org/3/library/pathlib.html](https://docs.python.org/3/library/pathlib.html) - Path operations
- Python argparse docs: [https://docs.python.org/3/library/argparse.html](https://docs.python.org/3/library/argparse.html) - CLI subcommands
- Current codebase: `/home/chris/develop/qm-nmr-calc/src/qm_nmr_calc/nwchem/geometry.py` - Existing RDKit patterns

### Secondary (MEDIUM confidence)
- NMReData initiative: [https://www.nmredata.org/](https://www.nmredata.org/) - SDF tag format for NMR data
- tqdm overview (2026): [https://www.analyticsvidhya.com/blog/2021/05/how-to-use-progress-bars-in-python/](https://www.analyticsvidhya.com/blog/2021/05/how-to-use-progress-bars-in-python/) - Progress bar best practices
- JSON Lines spec: [https://jsonlines.org/](https://jsonlines.org/) - JSONL format for streaming results
- Pydantic serialization: [https://docs.pydantic.dev/latest/concepts/serialization/](https://docs.pydantic.dev/latest/concepts/serialization/) - Model JSON handling
- pandas to_csv: [https://www.geeksforgeeks.org/python/how-to-append-pandas-dataframe-to-existing-csv-file/](https://www.geeksforgeeks.org/python/how-to-append-pandas-dataframe-to-existing-csv-file/) - CSV append mode

### Tertiary (LOW confidence)
- WebSearch-derived CLI patterns: [https://medium.com/@selvakumar-arumugapandian/command-line-subcommands-with-pythons-argparse-4dbac80f7110](https://medium.com/@selvakumar-arumugapandian/command-line-subcommands-with-pythons-argparse-4dbac80f7110) - Argparse subcommand examples (pattern verified with official docs)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - RDKit/Pydantic/argparse all verified in current dependencies or stdlib
- Architecture: HIGH - Patterns verified in current codebase (geometry.py, models.py) and official docs
- DELTA50 data: HIGH - Publication and SI confirmed available, format verified in PMC article
- NMReData format: MEDIUM - Standard exists but not required for internal benchmarking
- WP04 scaling factors: LOW - Not confirmed in current codebase, may need extraction from DELTA50 SI

**Research date:** 2026-01-21
**Valid until:** ~30 days (DELTA50 dataset stable, Python ecosystem stable)
