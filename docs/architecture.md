# Technical Architecture

> This guide will be completed in Phase 28.

## Topics to Be Covered

This guide will document the system architecture for developers and contributors:

- **Full Stack Overview**
  - FastAPI web framework
  - Huey task queue with SQLite backend
  - NWChem quantum chemistry engine
  - RDKit cheminformatics toolkit
  - 3Dmol.js molecule visualization

- **Data Flow Diagrams**
  - Request handling pipeline
  - Job submission to completion flow
  - Result generation and delivery

- **Job Lifecycle States**
  - State machine definition
  - Transition triggers
  - Error and recovery states

- **File Storage Structure**
  - Job directory layout
  - Scratch file management
  - Output artifact organization

- **Conformer Ensemble Pipeline**
  - CREST integration
  - Conformer selection and filtering
  - Parallel DFT execution
  - Boltzmann aggregation

- **CSS Architecture**
  - Layer organization
  - Component library
  - Design tokens and theming

---

For a high-level overview, see the [Architecture section](../README.md#architecture-overview) in the main README.
