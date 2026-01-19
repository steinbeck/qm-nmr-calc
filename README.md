# QM NMR Calculator

Asynchronous web service for NMR quantum mechanical calculations using ISiCLE and NWChem.

## Installation

```bash
uv sync
```

## Usage

Start the Huey consumer for background processing:

```bash
uv run huey_consumer.py qm_nmr_calc.queue.huey
```

Start the API server:

```bash
uv run uvicorn qm_nmr_calc.api.app:app
```

## Development

See `.planning/` for project planning and documentation.
