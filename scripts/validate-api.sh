#!/bin/bash
# Validation script for qm-nmr-calc API container
# Run inside container: docker run --rm qm-nmr-calc-api:test /app/scripts/validate-api.sh

set -e  # Exit on first error

echo "=== API Container Validation ==="
echo ""

# 1. Python environment validation
echo "--- Testing Python Environment ---"
which python || { echo "FAIL: python not in PATH"; exit 1; }
echo "Python found at: $(which python)"
python --version

# 2. Uvicorn validation
echo ""
echo "--- Testing Uvicorn ---"
which uvicorn || { echo "FAIL: uvicorn not in PATH"; exit 1; }
echo "Uvicorn found at: $(which uvicorn)"
uvicorn --version

# 3. FastAPI app import validation
echo ""
echo "--- Testing FastAPI App Import ---"
python -c "from qm_nmr_calc.api.app import app; print('App title:', app.title)" || {
    echo "FAIL: Cannot import FastAPI app"
    exit 1
}
echo "PASS: FastAPI app imports successfully"

# 4. Static files validation
echo ""
echo "--- Testing Static Files ---"
STATIC_DIR="/app/src/qm_nmr_calc/api/static"
if [ -d "$STATIC_DIR" ]; then
    CSS_COUNT=$(find "$STATIC_DIR" -name "*.css" | wc -l)
    echo "Static directory exists: $STATIC_DIR"
    echo "CSS files found: $CSS_COUNT"
    if [ "$CSS_COUNT" -gt 0 ]; then
        echo "PASS: Static files present"
    else
        echo "FAIL: No CSS files found"
        exit 1
    fi
else
    echo "FAIL: Static directory not found: $STATIC_DIR"
    exit 1
fi

# 5. Templates validation
echo ""
echo "--- Testing Templates ---"
TEMPLATE_DIR="/app/src/qm_nmr_calc/api/templates"
if [ -d "$TEMPLATE_DIR" ]; then
    HTML_COUNT=$(find "$TEMPLATE_DIR" -name "*.html" | wc -l)
    echo "Templates directory exists: $TEMPLATE_DIR"
    echo "HTML files found: $HTML_COUNT"
    if [ "$HTML_COUNT" -gt 0 ]; then
        echo "PASS: Templates present"
    else
        echo "FAIL: No HTML templates found"
        exit 1
    fi
else
    echo "FAIL: Templates directory not found: $TEMPLATE_DIR"
    exit 1
fi

# 6. Data directory validation
echo ""
echo "--- Testing Data Directory ---"
DATA_DIR="/app/data"
if [ -d "$DATA_DIR" ]; then
    if [ -w "$DATA_DIR" ]; then
        echo "Data directory exists and writable: $DATA_DIR"
        echo "PASS: Data directory accessible"
    else
        echo "FAIL: Data directory not writable: $DATA_DIR"
        exit 1
    fi
else
    echo "FAIL: Data directory not found: $DATA_DIR"
    exit 1
fi

# 7. Non-root user validation
echo ""
echo "--- Testing User ---"
CURRENT_USER=$(whoami)
CURRENT_UID=$(id -u)
echo "Running as: $CURRENT_USER (UID: $CURRENT_UID)"
if [ "$CURRENT_UID" -ne 0 ]; then
    echo "PASS: Running as non-root user"
else
    echo "FAIL: Running as root (UID 0)"
    exit 1
fi

# 8. Health endpoint validation (start uvicorn in background)
echo ""
echo "--- Testing Health Endpoint ---"
uvicorn qm_nmr_calc.api.app:app --host 0.0.0.0 --port 8000 &
UVICORN_PID=$!
sleep 3  # Wait for startup

# Test liveness endpoint
HTTP_CODE=$(curl -s -o /dev/null -w "%{http_code}" http://localhost:8000/health || echo "000")
if [ "$HTTP_CODE" = "200" ]; then
    echo "PASS: /health returns 200"
else
    echo "FAIL: /health returned $HTTP_CODE (expected 200)"
    kill $UVICORN_PID 2>/dev/null
    exit 1
fi

# Test static file access
HTTP_CODE=$(curl -s -o /dev/null -w "%{http_code}" http://localhost:8000/static/css/base.css || echo "000")
if [ "$HTTP_CODE" = "200" ]; then
    echo "PASS: Static CSS accessible"
else
    echo "FAIL: Static CSS returned $HTTP_CODE (expected 200)"
    kill $UVICORN_PID 2>/dev/null
    exit 1
fi

# Cleanup
kill $UVICORN_PID 2>/dev/null

echo ""
echo "=== All Validations Passed ==="
