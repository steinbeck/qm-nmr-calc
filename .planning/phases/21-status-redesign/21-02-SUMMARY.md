# Summary: 21-02 Status Page Template Update

## Overview

Updated the status.html template to use the new bento grid layout with glass cards, integrated the step tracker component, and updated JavaScript to build the step tracker dynamically. Includes additional fixes made during v2.0.1 testing.

## Deliverables

| Artifact | Location | Purpose |
|----------|----------|---------|
| status.html | src/qm_nmr_calc/api/templates/status.html | Status page with bento grid layout |

## Commits

| Task | Commit | Files |
|------|--------|-------|
| Task 1: Update template structure | b8139f8 | status.html |
| Task 2: Fix step tracker IDs | 1f26342 | status.html |
| Post-plan fixes: Step tracker improvements | 07f6261 | status.html |

## Key Decisions

1. **Step IDs match backend**: Step tracker uses backend step IDs (generating_conformers, optimizing_conformers, calculating_nmr, averaging_shifts, post_processing) for accurate progress display
2. **Sum durations for repeated steps**: JavaScript sums duration_seconds for steps that appear multiple times (e.g., optimizing_conformers runs once per conformer)
3. **Removed conformer details dropdown**: The accordion table showing per-conformer status was removed as it only updated on job completion (when page redirects anyway)

## Components Used

- **Bento Grid**: 6-column grid with span-3 and span-6 items
- **Glass Cards**: Job details, calculation progress, conformer progress, error display
- **Viewer Card (compact)**: 3D molecule preview at 250px height
- **Step Tracker**: Horizontal timeline with complete/active/pending states
- **Status Badge**: Header badge with data-status attribute for styling
- **Conformer Progress**: Native progress bar with ETA display (ensemble only)
- **Error Card**: Red variant with exclamation icon

## Verification

```bash
grep -c "bento-grid" src/qm_nmr_calc/api/templates/status.html    # 6
grep -c "step-tracker" src/qm_nmr_calc/api/templates/status.html  # 10
grep -c "glass-card" src/qm_nmr_calc/api/templates/status.html    # 16
grep -c "viewer-card" src/qm_nmr_calc/api/templates/status.html   # 3
```

## Human Verification

Verified during v2.0.1 hexanol calculation:
- Bento grid layout displays correctly
- Step tracker shows accurate progress with timing
- 3D molecule preview renders
- Conformer progress bar updates for ensemble jobs
- Status badge shows running/complete states

## Next

Phase 21 complete. Proceed to Phase 22 (Responsive and Layout Polish).
