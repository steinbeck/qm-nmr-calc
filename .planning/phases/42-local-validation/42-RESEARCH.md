# Phase 42: Local Validation - Research

**Researched:** 2026-02-04
**Domain:** ARM64 container validation for computational chemistry (NWChem, xTB, CREST)
**Confidence:** HIGH

## Summary

This phase validates the ARM64 worker container created in Phase 41 by running actual computational chemistry calculations on Apple Silicon. The validation must go beyond the simple binary/import checks in the existing `validate-worker-arm64.sh` script to actually execute calculations and verify numerical output.

The key challenge is establishing appropriate numerical tolerances for comparing ARM64 results to x86. Different BLAS implementations (OpenBLAS on ARM64 vs potentially MKL-optimized on x86) can produce slightly different floating-point results due to different algorithm implementations and instruction ordering, but these differences should be at machine-precision level (1e-10 or smaller for total energies). For chemical shifts, tolerances of 0.1-0.5 ppm for 1H and 1-2 ppm for 13C are appropriate given inherent DFT method uncertainty.

**Primary recommendation:** Create a staged validation script that builds the ARM64 image locally, runs it with `docker run --rm`, and executes progressively more complex calculations: (1) binary validation, (2) simple NWChem DFT optimization, (3) NWChem NMR shielding, (4) xTB energy, (5) CREST conformer search, (6) full pipeline comparison.

## Standard Stack

### Core (Validation Tools)

| Tool | Version | Purpose | Notes |
|------|---------|---------|-------|
| Docker Desktop | Latest | Build/run ARM64 images | Native ARM64 on Apple Silicon |
| bash | Built-in | Validation script | Portable across systems |
| pytest | >=9.0.2 | Test runner | Already in dev dependencies |

### Supporting

| Tool | Purpose | When to Use |
|------|---------|-------------|
| `docker buildx` | Build ARM64 images | Building on native Apple Silicon |
| `docker run --rm` | Ephemeral container execution | Running validation tests |
| `timeout` command | Prevent hung tests | Wrap long-running calculations |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| bash validation script | pytest-docker plugin | Simpler approach without extra dependency |
| Manual comparison | Reference output files | More robust but requires maintaining reference data |
| Single validation script | Multiple focused scripts | Single script easier for user execution |

## Architecture Patterns

### Validation Script Structure

```bash
#!/bin/bash
# scripts/validate-worker-arm64-full.sh
# Full computational validation (not just binary checks)

set -e

IMAGE_NAME="${1:-qm-nmr-calc-worker-arm64:test}"
TIMEOUT_NWCHEM=300  # 5 min for DFT
TIMEOUT_CREST=600   # 10 min for CREST

# 1. Build image (if needed)
# 2. Run binary validation (existing script)
# 3. NWChem geometry optimization test
# 4. NWChem NMR shielding test
# 5. xTB energy calculation test
# 6. CREST conformer search test
# 7. Python application test
# 8. Full pipeline test
```

### Test Data Pattern

```
tests/
  fixtures/
    arm64_validation/
      methane.xyz           # Simple test molecule
      ethanol_optimized.xyz # Pre-optimized geometry for NMR
      expected_nwchem.out   # Reference output patterns
```

### Docker Run Pattern for Tests

```bash
# Run test inside container with timeout
docker run --rm \
  --platform linux/arm64 \
  -v "$(pwd)/test_data:/app/test_data:ro" \
  "$IMAGE_NAME" \
  timeout 300 bash -c "cd /tmp && nwchem /app/test_data/test.nw"
```

### Anti-Patterns to Avoid

- **Testing without timeouts:** DFT calculations can hang; always use timeout
- **Testing only binaries:** `which nwchem` passing doesn't mean calculations work
- **Ignoring exit codes:** NWChem returns 0 even on some errors; check output
- **Hardcoded paths:** Use environment variables for portability
- **Running tests as root:** Match container's non-root user permissions

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| NMR validation | Custom parsing | Existing `parse_shielding_output()` | Already handles NWChem format |
| CREST validation | Custom parsing | Existing `parse_crest_ensemble()` | Already handles multi-XYZ format |
| Energy comparison | String matching | pytest.approx() | Handles floating-point tolerance |
| Timeout handling | Signal handlers | `timeout` command or subprocess timeout | Battle-tested, cross-platform |

**Key insight:** The codebase already has parsers for all output formats. Validation should reuse these parsers and check that outputs match expected patterns.

## Common Pitfalls

### Pitfall 1: SIGILL on ARM64

**What goes wrong:** NWChem crashes with "Illegal instruction" (SIGILL)
**Why it happens:** Binary compiled with x86-specific instructions, or QEMU emulation issues
**How to avoid:** Use conda-forge packages (already done in Phase 41), run on native ARM64
**Warning signs:** Immediate crash on `nwchem` execution, segfault during BLAS operations
**Validation:** First test should be simplest possible calculation (single-point energy)

### Pitfall 2: Basis Set Not Found

**What goes wrong:** NWChem errors: "Could not find basis set library"
**Why it happens:** NWCHEM_BASIS_LIBRARY not set or pointing to wrong location
**How to avoid:** Verify conda activation scripts set the variable correctly
**Warning signs:** "basis <name> does not exist" in output
**Validation:** Check `$NWCHEM_BASIS_LIBRARY` before running calculations

### Pitfall 3: OpenMPI Issues in Container

**What goes wrong:** NWChem hangs or crashes with MPI errors
**Why it happens:** MPI inside container needs special handling, shm_size limits
**How to avoid:** Use `--shm-size=512m` in docker run, consider single-process for validation
**Warning signs:** "MPI initialization failed", hangs at "Initializing Global Operations"
**Validation:** Start with `NWCHEM_NPROC=1` to isolate issues

### Pitfall 4: Stack Overflow on CREST/xTB

**What goes wrong:** CREST crashes on molecules with 30+ atoms
**Why it happens:** Default stack size insufficient, OMP_STACKSIZE not set
**How to avoid:** Verify OMP_STACKSIZE=2G is set in container environment
**Warning signs:** SIGSEGV during conformer search, works on small molecules
**Validation:** Test CREST on ethanol (9 atoms) then larger molecule

### Pitfall 5: Numerical Tolerance Too Strict

**What goes wrong:** Tests fail due to floating-point differences between architectures
**Why it happens:** Different BLAS implementations have different rounding behavior
**How to avoid:** Use appropriate tolerances (see table below)
**Warning signs:** Failures with differences at 1e-6 or smaller
**Validation:** Check that differences are consistent across multiple runs

### Pitfall 6: X11/Drawing Errors

**What goes wrong:** RDKit mol drawing fails in headless container
**Why it happens:** X11 libraries missing or DISPLAY not set
**How to avoid:** Verify xorg-libxrender, xorg-libxext installed (Phase 41 did this)
**Warning signs:** "cannot open display" errors
**Validation:** Run `python -c "from rdkit.Chem.Draw import MolToImage; ..."`

## Code Examples

### Minimal NWChem Test Input

```
# test_methane_opt.nw
# Minimal DFT optimization for validation
title "Methane ARM64 Validation"
memory 500 mb
geometry units angstrom
  C    0.000000    0.000000    0.000000
  H    0.629118    0.629118    0.629118
  H   -0.629118   -0.629118    0.629118
  H   -0.629118    0.629118   -0.629118
  H    0.629118   -0.629118   -0.629118
end
basis
  * library 6-31G*
end
dft
  xc b3lyp
  iterations 100
end
task dft optimize
```

### NWChem NMR Shielding Test

```
# test_methane_nmr.nw
title "Methane NMR ARM64 Validation"
memory 500 mb
geometry units angstrom
  # Pre-optimized geometry
  C    0.000000    0.000000    0.000000
  H    0.629118    0.629118    0.629118
  H   -0.629118   -0.629118    0.629118
  H   -0.629118    0.629118   -0.629118
  H    0.629118   -0.629118   -0.629118
end
basis
  * library 6-311+G(2d,p)
end
dft
  xc b3lyp
  iterations 100
end
property
  shielding
end
task dft property
```

### xTB Energy Test

```bash
# Inside container
cd /tmp
cat > test.xyz << EOF
5

C    0.000000    0.000000    0.000000
H    0.629118    0.629118    0.629118
H   -0.629118   -0.629118    0.629118
H   -0.629118    0.629118   -0.629118
H    0.629118   -0.629118   -0.629118
EOF

xtb test.xyz --gfn2 --verbose
# Check exit code and parse energy from output
```

### CREST Conformer Test

```bash
# Inside container
cd /tmp
cat > ethanol.xyz << EOF
9

C   -0.000  0.000  0.000
C    1.520  0.000  0.000
O   -0.537  1.310  0.000
H   -0.363 -0.513  0.889
H   -0.363 -0.513 -0.889
H    1.893  0.513 -0.889
H    1.893  0.513  0.889
H    1.893 -1.026  0.000
H   -0.100  1.820  0.720
EOF

crest ethanol.xyz --gfn2 --ewin 6.0 -T 2 --noreftopo
# Check crest_conformers.xyz exists and has multiple conformers
```

### Python Validation Test

```python
# test_arm64_validation.py
"""Integration tests for ARM64 container validation."""
import subprocess
import tempfile
from pathlib import Path
import pytest

# Skip if not running inside ARM64 container
pytestmark = pytest.mark.skipif(
    # Check conditions
    reason="ARM64 container tests"
)

class TestNWChemARM64:
    """NWChem validation on ARM64."""

    def test_geometry_optimization_completes(self, tmp_path):
        """DFT geometry optimization completes without SIGILL."""
        input_text = """..."""  # Methane input
        input_file = tmp_path / "test.nw"
        input_file.write_text(input_text)

        result = subprocess.run(
            ["nwchem", str(input_file)],
            capture_output=True,
            timeout=300,
            text=True
        )

        assert result.returncode == 0, f"NWChem failed: {result.stderr}"
        assert "optimization converged" in result.stdout.lower()

    def test_nmr_shielding_produces_valid_output(self, tmp_path):
        """NMR shielding calculation produces parseable output."""
        # Similar pattern...
```

## Numerical Tolerance Guidelines

| Quantity | Tolerance | Rationale |
|----------|-----------|-----------|
| Total energy (Hartree) | 1e-6 | DFT convergence threshold |
| Relative energy (kcal/mol) | 0.01 | Derived from Hartree tolerance |
| Geometry (Angstrom) | 1e-4 | Optimization convergence |
| Shielding constant (ppm) | 0.1 | Below method accuracy |
| 1H chemical shift | 0.5 ppm | Typical DFT accuracy ~0.2-0.5 ppm |
| 13C chemical shift | 2.0 ppm | Typical DFT accuracy ~2-5 ppm |
| CREST conformer count | Exact | Should be reproducible |
| Boltzmann weights | 0.01 | Energy-derived, sensitive |

**Key insight:** ARM64 vs x86 differences should be well below method accuracy. If differences exceed these tolerances, it indicates a real problem (wrong basis set, failed convergence, etc.), not architecture effects.

## State of the Art

| Old Approach | Current Approach | Why Different |
|--------------|------------------|---------------|
| Manual testing | Scripted validation | Reproducibility, CI integration |
| Binary-only checks | Full calculation tests | Binaries can exist but fail on execution |
| Single test molecule | Progressive complexity | Isolate failure modes |
| No timeout | Timeout on all calculations | Prevent hung tests |
| x86 reference values | Tolerance-based comparison | Architecture-independent |

## Open Questions

1. **Exact Reference Values**
   - What we know: Tolerances defined above should be sufficient
   - What's unclear: Whether to store x86 reference values or just validate output format
   - Recommendation: Start with format validation; add reference values if needed

2. **Large Molecule Performance**
   - What we know: CREST scales with molecule size, needs OMP_STACKSIZE
   - What's unclear: Largest molecule that reliably completes in reasonable time
   - Recommendation: Test with ethanol (9 atoms), limit validation to <20 atoms

3. **Multi-Process NWChem**
   - What we know: Phase 41 sets NWCHEM_NPROC=4
   - What's unclear: Whether MPI works correctly in ARM64 container
   - Recommendation: Test with NPROC=1 first, then verify NPROC>1

## Success Criteria Verification

| Criterion | How to Verify |
|-----------|---------------|
| NWChem DFT geometry optimization completes without SIGILL or crashes | Run optimization, check exit code 0, parse "optimization converged" |
| NWChem NMR shielding calculation produces valid output | Run property task, parse with existing `parse_shielding_output()`, verify H/C values |
| xTB energy calculation produces expected format output | Run `xtb --gfn2`, check exit code, verify energy line in output |
| CREST conformer search generates multiple conformers | Run `crest` on ethanol, verify `crest_conformers.xyz` has >1 structure |
| Python application imports and task queue initializes | `python -c "from qm_nmr_calc.queue import huey"` |
| Full NMR prediction pipeline produces results matching x86 within tolerance | Run end-to-end on test molecule, compare shifts within tolerance |

## Sources

### Primary (HIGH confidence)
- [NWChem Geometry Optimization](https://nwchemgit.github.io/Geometry-Optimization.html) - Convergence parameters
- [NWChem DFT](https://nwchemgit.github.io/Density-Functional-Theory-for-Molecules.html) - DFT module documentation
- [CREST Conformational Sampling](https://crest-lab.github.io/crest-docs/page/examples/example_1.html) - Output format, expected files
- [Docker Multi-Platform](https://docs.docker.com/build/building/multi-platform/) - `--platform` flag usage
- Phase 41 RESEARCH.md and SUMMARY.md - ARM64 container configuration

### Secondary (MEDIUM confidence)
- [OpenBLAS ARM64 Issues](https://github.com/OpenMathLib/OpenBLAS/issues/2814) - Tuning for Apple chips
- [Psi4 Apple Silicon Notes](https://github.com/psi4/psi4/issues/2333) - Computational chemistry on ARM64
- [CREST GitHub](https://github.com/crest-lab/crest) - Output file patterns

### Codebase (HIGH confidence - verified)
- `tests/test_nwchem_integration.py` - Existing NWChem test patterns
- `tests/test_crest_generator.py` - CREST parsing tests
- `tests/test_xtb_ranking.py` - xTB validation tests
- `scripts/validate-worker-arm64.sh` - Existing binary validation (to extend)
- `src/qm_nmr_calc/nwchem/output_parser.py` - Shielding parsing

## Metadata

**Confidence breakdown:**
- Validation approach: HIGH - Based on existing test patterns and documented tools
- Numerical tolerances: MEDIUM - Extrapolated from DFT accuracy, may need tuning
- Docker patterns: HIGH - Based on official Docker documentation

**Research date:** 2026-02-04
**Valid until:** 2026-03-04 (30 days - stable domain)
