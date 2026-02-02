#!/bin/bash
# Full NWChem DFT test inside container
# This verifies NWChem can actually run calculations, not just respond to --version
# CRITICAL: Run with --shm-size=512m to provide sufficient shared memory for MPI

set -e

echo "=== Full NWChem DFT Test ==="

# Create a simple test input for water molecule optimization
WORKDIR=$(mktemp -d)
cd "$WORKDIR"

cat > water.nw << 'EOF'
start water
title "Water geometry optimization"
memory total 1 gb

geometry units angstrom
  O  0.000  0.000  0.117
  H  0.000  0.756 -0.469
  H  0.000 -0.756 -0.469
end

basis
  * library 6-31g*
end

driver
  maxiter 50
end

task dft optimize
EOF

echo "Running NWChem optimization on water molecule..."
echo "Input file: water.nw"
echo ""

# Run with MPI (single process for test, using --bind-to none for container compatibility)
mpirun --bind-to none -n 1 nwchem water.nw > water.out 2> water.err || {
    echo "FAIL: NWChem calculation failed"
    echo "=== STDERR ==="
    cat water.err
    echo "=== STDOUT (last 50 lines) ==="
    tail -50 water.out
    exit 1
}

# Check for successful completion
if grep -q "Total DFT energy" water.out; then
    echo "PASS: NWChem DFT calculation completed"
    grep "Total DFT energy" water.out
else
    echo "FAIL: DFT energy not found in output"
    tail -50 water.out
    exit 1
fi

# Check geometry optimization converged
if grep -q "Optimization converged" water.out; then
    echo "PASS: Geometry optimization converged"
else
    echo "WARNING: Optimization may not have converged (check output)"
fi

# Cleanup
cd /
rm -rf "$WORKDIR"

echo ""
echo "=== NWChem DFT Test Passed ==="
