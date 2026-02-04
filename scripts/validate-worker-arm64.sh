#!/bin/bash
# Validation script for qm-nmr-calc ARM64 worker container
# Run inside container: docker run --rm qm-nmr-calc-worker-arm64:test /app/scripts/validate-worker-arm64.sh

set -e  # Exit on first error

echo "=== ARM64 Worker Container Validation ==="
echo ""

# 0. Architecture check (warning, not failure)
echo "--- Checking Architecture ---"
ARCH=$(uname -m)
echo "Architecture: $ARCH"
if [ "$ARCH" != "aarch64" ]; then
    echo "WARNING: Not running on ARM64 (aarch64). Detected: $ARCH"
    echo "         This may be QEMU emulation or wrong architecture."
    echo "         Continuing with validation..."
else
    echo "PASS: Running on ARM64 (aarch64)"
fi
echo ""

# 1. NWChem validation
echo "--- Testing NWChem ---"
which nwchem || { echo "FAIL: nwchem not in PATH"; exit 1; }
echo "NWChem found at: $(which nwchem)"

# NWChem doesn't have --version flag; verify by checking it's executable
if [ -x "$(which nwchem)" ]; then
    echo "PASS: NWChem binary exists and is executable"
else
    echo "FAIL: NWChem not executable"
    exit 1
fi

# Check NWCHEM_BASIS_LIBRARY is set (by conda activation)
if [ -n "$NWCHEM_BASIS_LIBRARY" ]; then
    echo "NWCHEM_BASIS_LIBRARY: $NWCHEM_BASIS_LIBRARY"
    if [ -d "$NWCHEM_BASIS_LIBRARY" ]; then
        BASIS_COUNT=$(ls -1 "$NWCHEM_BASIS_LIBRARY" 2>/dev/null | wc -l)
        echo "Basis sets found: $BASIS_COUNT files"
        echo "PASS: Basis library exists with files"
    else
        echo "FAIL: NWCHEM_BASIS_LIBRARY directory does not exist"
        exit 1
    fi
else
    echo "FAIL: NWCHEM_BASIS_LIBRARY not set (conda activation may have failed)"
    exit 1
fi
echo ""

# 2. xTB validation
echo "--- Testing xTB ---"
which xtb || { echo "FAIL: xtb not in PATH"; exit 1; }
echo "xTB found at: $(which xtb)"
xtb --version 2>&1 | head -3 || { echo "FAIL: xtb --version failed"; exit 1; }
echo "PASS: xTB responds"
echo ""

# 3. CREST validation
echo "--- Testing CREST ---"
which crest || { echo "FAIL: crest not in PATH"; exit 1; }
echo "CREST found at: $(which crest)"
crest --version 2>&1 | head -3 || { echo "FAIL: crest --version failed"; exit 1; }
echo "PASS: CREST responds"
echo ""

# 4. Python environment validation
echo "--- Testing Python Environment ---"
which python || { echo "FAIL: python not in PATH"; exit 1; }
echo "Python found at: $(which python)"
python --version
PYTHON_VERSION=$(python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
if [ "$PYTHON_VERSION" = "3.11" ]; then
    echo "PASS: Python 3.11"
else
    echo "WARNING: Python version is $PYTHON_VERSION (expected 3.11)"
fi
echo ""

# 5. RDKit validation (from conda-forge)
echo "--- Testing RDKit ---"
python -c "
from rdkit import Chem
from rdkit import rdBase
mol = Chem.MolFromSmiles('CCO')
print(f'RDKit version: {rdBase.rdkitVersion}')
print(f'Test molecule (ethanol): {mol.GetNumAtoms()} atoms')
" || {
    echo "FAIL: RDKit not working"
    exit 1
}
echo "PASS: RDKit works"
echo ""

# 6. Huey import validation
echo "--- Testing Huey Import ---"
python -c "from qm_nmr_calc.queue import huey; print('Huey instance:', huey.name)" || {
    echo "FAIL: Cannot import Huey from qm_nmr_calc.queue"
    exit 1
}
echo "PASS: Huey imports successfully"
echo ""

# 7. Task import validation
echo "--- Testing Task Import ---"
python -c "from qm_nmr_calc.tasks import run_nmr_task; print('Task:', run_nmr_task)" || {
    echo "FAIL: Cannot import tasks from qm_nmr_calc.tasks"
    exit 1
}
echo "PASS: Tasks import successfully"
echo ""

# 8. Environment variables check
echo "--- Checking Environment Variables ---"
echo "OMP_STACKSIZE: ${OMP_STACKSIZE:-NOT SET}"
echo "OMP_NUM_THREADS: ${OMP_NUM_THREADS:-NOT SET}"
echo "OPENBLAS_NUM_THREADS: ${OPENBLAS_NUM_THREADS:-NOT SET}"
echo "NWCHEM_NPROC: ${NWCHEM_NPROC:-NOT SET}"

if [ -n "$OMP_STACKSIZE" ] && [ -n "$OMP_NUM_THREADS" ]; then
    echo "PASS: OpenMP variables set"
else
    echo "WARNING: Some OpenMP variables not set"
fi
echo ""

echo "=== All Validations Passed ==="
exit 0
