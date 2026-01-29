# Plan 17-05 Summary: Ensemble Results Display

## What Was Built

Added ensemble metadata display and clarified shift averaging on results page.

### Ensemble Metadata Section
- Shows "Ensemble Averaging Summary" for ensemble jobs
- Displays: conformers used, total generated, method, temperature, energy range
- Top contributing conformers table with ID, population, relative energy
- Clear note: "The chemical shifts displayed below are Boltzmann-weighted averages"
- Warning banner for CREST fallback scenarios

### Shift Data Table
- Expandable "View Shift Data Table" section
- Shows averaged 1H and 13C shifts in tabular format
- Column headers clarify "Averaged Shift"

### 3D Viewer Clarifications
- Legend updated to "1H averaged shifts (ppm)" and "13C averaged shifts (ppm)"
- Conformer info changed to "Geometry view — this conformer has X% population"
- Added note: "Shift labels show Boltzmann-averaged values. Use dropdown to view conformer geometries."

## Files Modified

- `src/qm_nmr_calc/api/templates/results.html` - Ensemble metadata, shift tables, viewer legend

## Key UX Decisions

- Shift labels on 3D viewer always show averaged values (not per-conformer)
- Conformer selector changes geometry view only, labels remain averaged
- Clear visual separation between geometry view and shift data

## Verification

- Ensemble jobs show metadata section with conformer breakdown
- Single-conformer jobs do NOT show ensemble sections
- Legend and notes clearly explain what the shift values represent
- End-to-end test with real ensemble job confirmed working

## Duration

~20 minutes (including iterative UX refinement)
