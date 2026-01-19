#!/usr/bin/env python
"""Run the QM NMR Calculator API server."""
import sys
from pathlib import Path

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import uvicorn


def main():
    """Run the API server with auto-reload for development."""
    uvicorn.run(
        "qm_nmr_calc.api.app:app",
        host="0.0.0.0",
        port=8000,
        reload=True,
        reload_dirs=["src"],
    )


if __name__ == "__main__":
    main()
