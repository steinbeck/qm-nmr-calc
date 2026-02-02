#!/bin/bash
# Validation script for qm-nmr-calc worker container
# Run inside container: docker run --rm qm-nmr-calc-worker:test /app/scripts/validate-worker.sh

set -e  # Exit on first error

echo "=== Worker Container Validation ==="
echo ""

# 1. NWChem validation
echo "--- Testing NWChem ---"
which nwchem || { echo "FAIL: nwchem not in PATH"; exit 1; }
echo "NWChem found at: $(which nwchem)"

# NWChem doesn't have --version flag; verify by checking it's executable and basis library exists
if [ -x "$(which nwchem)" ] && [ -d "${NWCHEM_BASIS_LIBRARY:-/opt/nwchem/src/basis/libraries}" ]; then
    echo "PASS: NWChem binary exists and is executable"
    echo "Basis library: ${NWCHEM_BASIS_LIBRARY:-/opt/nwchem/src/basis/libraries}"
else
    echo "FAIL: NWChem not properly configured"
    exit 1
fi

# 2. xTB validation
echo ""
echo "--- Testing xTB ---"
which xtb || { echo "FAIL: xtb not in PATH"; exit 1; }
echo "xTB found at: $(which xtb)"
xtb --version || { echo "FAIL: xtb --version failed"; exit 1; }
echo "PASS: xTB responds"

# 3. CREST validation
echo ""
echo "--- Testing CREST ---"
which crest || { echo "FAIL: crest not in PATH"; exit 1; }
echo "CREST found at: $(which crest)"
crest --version || { echo "FAIL: crest --version failed"; exit 1; }
echo "PASS: CREST responds"

# 4. Python environment validation
echo ""
echo "--- Testing Python Environment ---"
which python || { echo "FAIL: python not in PATH"; exit 1; }
echo "Python found at: $(which python)"
python --version

# 5. Huey import validation
echo ""
echo "--- Testing Huey Import ---"
python -c "from qm_nmr_calc.queue import huey; print('Huey instance:', huey.name)" || {
    echo "FAIL: Cannot import Huey from qm_nmr_calc.queue"
    exit 1
}
echo "PASS: Huey imports successfully"

# 6. Task import validation
echo ""
echo "--- Testing Task Import ---"
python -c "from qm_nmr_calc.tasks import run_nmr_task; print('Task:', run_nmr_task)" || {
    echo "FAIL: Cannot import tasks from qm_nmr_calc.tasks"
    exit 1
}
echo "PASS: Tasks import successfully"

# 7. RDKit validation (needed for conformer generation)
echo ""
echo "--- Testing RDKit ---"
python -c "from rdkit import Chem; m = Chem.MolFromSmiles('CCO'); print('Atoms:', m.GetNumAtoms())" || {
    echo "FAIL: RDKit not working"
    exit 1
}
echo "PASS: RDKit works"

echo ""
echo "=== All Validations Passed ==="
