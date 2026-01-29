# Plan 17-03 Summary: Web Form and Status Page

## What Was Built

Extended web form submission and status page to support ensemble mode parameters.

### Web Form (submit.html)
- Added conformer mode dropdown (Ensemble/Single) with ensemble as default
- Added conformer method dropdown (RDKit KDG/CREST) with dynamic CREST availability
- Mode/method dropdowns always visible per CONTEXT.md decision
- Form preserves values on validation errors

### Status Page (status.html)
- Added conformer progress section showing X/N counts
- Displays current step with conformer count (e.g., "Optimizing Conformers (3/6)")
- Shows method used and any warnings (e.g., CREST fallback)

### Web Router (web.py)
- Updated submit_job to accept conformer_mode, conformer_method, max_conformers
- Dispatch to run_ensemble_nmr_task for ensemble mode
- Default conformer_mode="ensemble" per CONTEXT.md decision

## Files Modified

- `src/qm_nmr_calc/api/templates/submit.html` - Mode/method dropdowns
- `src/qm_nmr_calc/api/templates/status.html` - Conformer progress display
- `src/qm_nmr_calc/api/routers/web.py` - Form handling with conformer params

## Verification

- Web form shows mode/method dropdowns
- Ensemble jobs dispatch correctly to run_ensemble_nmr_task
- Status page shows conformer progress during processing
- Single-conformer mode still works as before

## Duration

~15 minutes (including manual verification)
