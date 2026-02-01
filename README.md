# QM NMR Calculator

![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)
![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)

Predict NMR chemical shifts for organic molecules using quantum chemistry. Submit a molecule structure (SMILES, MOL, or SDF), and get back predicted 1H and 13C chemical shifts with interactive spectrum visualizations and a 3D structure viewer.

Built for chemists and researchers who need NMR predictions without managing complex computational workflows.

## Features

- **NMR chemical shift predictions** using DFT calculations with linear scaling calibration
- **Conformer ensembles** with Boltzmann-weighted averaging (via optional CREST)
- **Multiple solvents** including CHCl3, DMSO, and gas phase
- **Web interface** for chemists - no command line required
- **REST API** for automated workflows and integrations
- **Interactive results** with spectrum plots and 3D viewer with shift labels
- **Draft and production modes** for quick checks vs publication-quality results

## Getting Started

See the **[Installation Guide](docs/installation.md)** for system dependencies, setup instructions, and a quick start walkthrough.

## Documentation

- **[NMR Methodology](docs/science.md)** - DP4+, linear scaling, and Boltzmann averaging
- **[Usage Guide](docs/usage.md)** - Web UI workflow and REST API reference
- **[Technical Architecture](docs/architecture.md)** - System design and data flow
- **[Library Documentation](docs/libraries.md)** - RDKit, NWChem, Huey integrations

## API Example

```bash
# Submit calculation
curl -X POST http://localhost:8000/api/v1/jobs \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO", "solvent": "chcl3"}'

# Check status
curl http://localhost:8000/api/v1/jobs/{job_id}

# Get results (chemical shifts JSON)
curl http://localhost:8000/api/v1/jobs/{job_id}/results

# Get spectrum plot
curl http://localhost:8000/api/v1/jobs/{job_id}/spectrum/1H.svg -o spectrum.svg
```

## Supported Solvents

| Code | Solvent | Use Case |
|------|---------|----------|
| `chcl3` | Chloroform (CDCl3) | Standard organic chemistry |
| `dmso` | DMSO (DMSO-d6) | Polar compounds |

## Calculation Presets

| Preset | NMR Basis | Use Case |
|--------|-----------|----------|
| `draft` | 6-31G* | Quick checks, debugging (~5 min) |
| `production` | 6-311+G(2d,p) | Publication quality (~30 min) |

## Development

```bash
# Run tests
uv run pytest tests/ -v

# Project structure
qm-nmr-calc/
├── src/qm_nmr_calc/
│   ├── api/          # FastAPI application
│   ├── nwchem/       # NWChem integration
│   └── ...           # Core modules
├── scripts/          # Startup scripts
├── tests/            # Test suite
└── docs/             # Documentation
```

See [docs/architecture.md](docs/architecture.md) for detailed system design.

## License

MIT

## Acknowledgments

- [NWChem](https://nwchemgit.github.io/) - Quantum chemistry package for DFT calculations
- [RDKit](https://www.rdkit.org/) - Cheminformatics toolkit for molecular handling
- [ISiCLE](https://github.com/pnnl/isicle) - Methodology inspiration from Pacific Northwest National Laboratory
- [3Dmol.js](https://3dmol.csb.pitt.edu/) - Interactive molecule visualization
- [Huey](https://huey.readthedocs.io/) - Task queue for async processing
