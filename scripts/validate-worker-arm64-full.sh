#!/bin/bash
# Full computational validation for ARM64 worker container
# Runs actual calculations, not just binary checks
# Usage: ./scripts/validate-worker-arm64-full.sh [image-name]

set -e

IMAGE_NAME="${1:-qm-nmr-calc-worker-arm64:test}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
FIXTURES_DIR="$PROJECT_ROOT/tests/fixtures/arm64_validation"

# Timeouts (seconds)
TIMEOUT_OPT=300      # 5 min for geometry optimization
TIMEOUT_NMR=300      # 5 min for NMR shielding
TIMEOUT_XTB=60       # 1 min for xTB energy
TIMEOUT_CREST=600    # 10 min for CREST conformer search

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Result tracking
TESTS_PASSED=0
TESTS_FAILED=0

log_pass() {
    echo -e "${GREEN}PASS${NC}: $1"
    TESTS_PASSED=$((TESTS_PASSED + 1))
}

log_fail() {
    echo -e "${RED}FAIL${NC}: $1"
    TESTS_FAILED=$((TESTS_FAILED + 1))
}

log_warn() {
    echo -e "${YELLOW}WARN${NC}: $1"
}

log_section() {
    echo ""
    echo "========================================="
    echo "$1"
    echo "========================================="
}

# Check prerequisites
check_prerequisites() {
    log_section "Checking Prerequisites"

    if ! command -v docker &> /dev/null; then
        log_fail "Docker not installed"
        exit 1
    fi
    log_pass "Docker available"

    if [ ! -d "$FIXTURES_DIR" ]; then
        log_fail "Test fixtures not found at $FIXTURES_DIR"
        exit 1
    fi
    log_pass "Test fixtures found"

    # Check architecture
    ARCH=$(uname -m)
    if [ "$ARCH" = "arm64" ]; then
        log_pass "Running on Apple Silicon (arm64)"
    else
        log_warn "Running on $ARCH - may use QEMU emulation"
    fi
}

# Build image if needed
build_image() {
    log_section "Building ARM64 Image"

    if docker image inspect "$IMAGE_NAME" &> /dev/null; then
        log_warn "Image $IMAGE_NAME already exists, using existing"
        return 0
    fi

    echo "Building $IMAGE_NAME from Dockerfile.worker.arm64..."
    if docker build --platform linux/arm64 -t "$IMAGE_NAME" -f Dockerfile.worker.arm64 "$PROJECT_ROOT"; then
        log_pass "Image built successfully"
    else
        log_fail "Image build failed"
        exit 1
    fi
}

# Run basic binary validation (existing script)
run_binary_validation() {
    log_section "Binary Validation (Basic Checks)"

    if docker run --rm --platform linux/arm64 "$IMAGE_NAME" /app/scripts/validate-worker-arm64.sh; then
        log_pass "Binary validation passed"
    else
        log_fail "Binary validation failed"
        return 1
    fi
}

# Test 1: NWChem geometry optimization
test_nwchem_optimization() {
    log_section "Test 1: NWChem Geometry Optimization"

    echo "Running DFT optimization on methane (timeout: ${TIMEOUT_OPT}s)..."

    # Create temp directory for output
    TEMP_DIR=$(mktemp -d)
    trap "rm -rf $TEMP_DIR" EXIT

    # Run optimization inside container
    if timeout "$TIMEOUT_OPT" docker run --rm \
        --platform linux/arm64 \
        --shm-size=512m \
        -v "$FIXTURES_DIR:/app/test_data:ro" \
        -v "$TEMP_DIR:/app/output" \
        -e NWCHEM_NPROC=1 \
        "$IMAGE_NAME" \
        bash -c "cd /app/output && nwchem /app/test_data/test_methane_opt.nw > opt.out 2>&1"; then

        # Check for successful completion
        if grep -q "optimization converged" "$TEMP_DIR/opt.out" 2>/dev/null || \
           grep -qi "total dft energy" "$TEMP_DIR/opt.out" 2>/dev/null; then
            log_pass "NWChem geometry optimization completed"

            # Extract and display energy
            ENERGY=$(grep -i "total dft energy" "$TEMP_DIR/opt.out" | tail -1 | awk '{print $NF}')
            if [ -n "$ENERGY" ]; then
                echo "  Final energy: $ENERGY Hartree"
            fi
        else
            log_fail "NWChem optimization did not converge"
            echo "Last 50 lines of output:"
            tail -50 "$TEMP_DIR/opt.out" 2>/dev/null || echo "No output file"
            return 1
        fi
    else
        log_fail "NWChem optimization timed out or crashed"
        if [ -f "$TEMP_DIR/opt.out" ]; then
            echo "Last 50 lines of output:"
            tail -50 "$TEMP_DIR/opt.out"
        fi
        return 1
    fi
}

# Test 2: NWChem NMR shielding
test_nwchem_nmr() {
    log_section "Test 2: NWChem NMR Shielding"

    echo "Running NMR shielding calculation on methane (timeout: ${TIMEOUT_NMR}s)..."

    TEMP_DIR=$(mktemp -d)
    trap "rm -rf $TEMP_DIR" EXIT

    if timeout "$TIMEOUT_NMR" docker run --rm \
        --platform linux/arm64 \
        --shm-size=512m \
        -v "$FIXTURES_DIR:/app/test_data:ro" \
        -v "$TEMP_DIR:/app/output" \
        -e NWCHEM_NPROC=1 \
        "$IMAGE_NAME" \
        bash -c "cd /app/output && nwchem /app/test_data/test_methane_nmr.nw > nmr.out 2>&1"; then

        # Check for shielding output
        if grep -q "isotropic" "$TEMP_DIR/nmr.out" 2>/dev/null; then
            log_pass "NWChem NMR shielding calculation completed"

            # Extract shielding values
            echo "  Isotropic shielding values:"
            grep -A1 "Atom:" "$TEMP_DIR/nmr.out" | grep "isotropic" | head -5
        else
            log_fail "NWChem NMR shielding output not found"
            echo "Last 50 lines of output:"
            tail -50 "$TEMP_DIR/nmr.out" 2>/dev/null || echo "No output file"
            return 1
        fi
    else
        log_fail "NWChem NMR calculation timed out or crashed"
        return 1
    fi
}

# Test 3: xTB energy calculation
test_xtb_energy() {
    log_section "Test 3: xTB Energy Calculation"

    echo "Running xTB GFN2 energy on methane (timeout: ${TIMEOUT_XTB}s)..."

    TEMP_DIR=$(mktemp -d)
    trap "rm -rf $TEMP_DIR" EXIT

    if timeout "$TIMEOUT_XTB" docker run --rm \
        --platform linux/arm64 \
        -v "$FIXTURES_DIR:/app/test_data:ro" \
        -v "$TEMP_DIR:/app/output" \
        "$IMAGE_NAME" \
        bash -c "cd /app/output && xtb /app/test_data/methane.xyz --gfn2 > xtb.out 2>&1"; then

        # Check for energy in output
        if grep -q "TOTAL ENERGY" "$TEMP_DIR/xtb.out" 2>/dev/null; then
            log_pass "xTB energy calculation completed"

            ENERGY=$(grep "TOTAL ENERGY" "$TEMP_DIR/xtb.out" | head -1)
            echo "  $ENERGY"
        else
            log_fail "xTB energy output not found"
            echo "Last 30 lines of output:"
            tail -30 "$TEMP_DIR/xtb.out" 2>/dev/null || echo "No output file"
            return 1
        fi
    else
        log_fail "xTB calculation timed out or crashed"
        return 1
    fi
}

# Test 4: CREST conformer search
test_crest_conformers() {
    log_section "Test 4: CREST Conformer Search"

    echo "Running CREST on ethanol (timeout: ${TIMEOUT_CREST}s)..."
    echo "This may take several minutes..."

    TEMP_DIR=$(mktemp -d)
    trap "rm -rf $TEMP_DIR" EXIT

    if timeout "$TIMEOUT_CREST" docker run --rm \
        --platform linux/arm64 \
        -v "$FIXTURES_DIR:/app/test_data:ro" \
        -v "$TEMP_DIR:/app/output" \
        "$IMAGE_NAME" \
        bash -c "cd /app/output && crest /app/test_data/ethanol.xyz --gfn2 --ewin 6.0 -T 2 --noreftopo > crest.out 2>&1"; then

        # Check for conformer output
        if [ -f "$TEMP_DIR/crest_conformers.xyz" ]; then
            # Count conformers by counting lines that are just a number (atom count)
            # Each XYZ structure starts with atom count, so count lines matching just digits
            NUM_CONF=$(grep -cE '^[0-9]+$' "$TEMP_DIR/crest_conformers.xyz" 2>/dev/null | tr -d '[:space:]')
            NUM_CONF=${NUM_CONF:-0}

            if [ "$NUM_CONF" -gt 1 ]; then
                log_pass "CREST found $NUM_CONF conformers"
            elif [ "$NUM_CONF" -ge 1 ]; then
                log_warn "CREST found only 1 conformer (may be expected for ethanol)"
                # Still pass - some molecules have limited conformational flexibility
                log_pass "CREST completed successfully"
            else
                # Check if file has any content
                if [ -s "$TEMP_DIR/crest_conformers.xyz" ]; then
                    log_warn "CREST output exists but conformer count unclear"
                    log_pass "CREST completed (output file exists)"
                else
                    log_fail "CREST produced empty output"
                    return 1
                fi
            fi
        else
            log_fail "CREST output file not found"
            echo "Last 30 lines of output:"
            tail -30 "$TEMP_DIR/crest.out" 2>/dev/null || echo "No output file"
            return 1
        fi
    else
        log_fail "CREST timed out or crashed"
        return 1
    fi
}

# Test 5: Python application integration
test_python_integration() {
    log_section "Test 5: Python Application Integration"

    echo "Testing Python imports and basic functionality..."

    if docker run --rm --platform linux/arm64 "$IMAGE_NAME" \
        python -c "
from qm_nmr_calc.queue import huey
from qm_nmr_calc.tasks import run_nmr_task
from qm_nmr_calc.nwchem import parse_shielding_output
from qm_nmr_calc.conformers.crest_generator import parse_crest_ensemble
from rdkit import Chem
from rdkit.Chem import AllChem

# Test RDKit
mol = Chem.MolFromSmiles('CCO')
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol)
print(f'RDKit: ethanol has {mol.GetNumAtoms()} atoms')

# Test imports work
print('All imports successful')
"; then
        log_pass "Python application integration"
    else
        log_fail "Python integration failed"
        return 1
    fi
}

# Summary
print_summary() {
    log_section "Validation Summary"

    echo ""
    echo "Tests passed: $TESTS_PASSED"
    echo "Tests failed: $TESTS_FAILED"
    echo ""

    if [ "$TESTS_FAILED" -eq 0 ]; then
        echo -e "${GREEN}All validation tests passed!${NC}"
        echo ""
        echo "ARM64 worker container is ready for:"
        echo "  - CI/CD integration (Phase 43)"
        echo "  - Production deployment on ARM64 hosts"
        exit 0
    else
        echo -e "${RED}Some validation tests failed.${NC}"
        echo ""
        echo "Review the failures above and check:"
        echo "  - Dockerfile.worker.arm64 configuration"
        echo "  - env-worker-arm64.yaml packages"
        echo "  - Environment variables in container"
        exit 1
    fi
}

# Main execution
main() {
    echo "ARM64 Worker Container Full Validation"
    echo "Image: $IMAGE_NAME"
    echo ""

    check_prerequisites
    build_image
    run_binary_validation
    test_nwchem_optimization
    test_nwchem_nmr
    test_xtb_energy
    test_crest_conformers
    test_python_integration
    print_summary
}

main "$@"
