---
phase: 04-results-delivery
plan: 02
subsystem: api
tags: [email, notifications, aiosmtplib, async, smtp]

# Dependency graph
requires:
  - phase: 03-nmr-calculations
    provides: Job completion signals, status storage
provides:
  - Opt-in email notification on job completion
  - Email notification on job failure
  - notifications module with async/sync wrappers
affects: [05-monitoring, 06-deployment]

# Tech tracking
tech-stack:
  added: [aiosmtplib, email-validator]
  patterns: [async email with sync wrapper for Huey, best-effort notification delivery]

key-files:
  created:
    - src/qm_nmr_calc/notifications.py
  modified:
    - pyproject.toml
    - src/qm_nmr_calc/models.py
    - src/qm_nmr_calc/api/schemas.py
    - src/qm_nmr_calc/api/routers/jobs.py
    - src/qm_nmr_calc/storage.py
    - src/qm_nmr_calc/queue.py

key-decisions:
  - "aiosmtplib for async SMTP, email-validator for Pydantic EmailStr"
  - "Best-effort email delivery (logs errors, never fails jobs)"
  - "Environment variables for all SMTP config (no hardcoded credentials)"
  - "Multipart email with plain text and HTML versions"

patterns-established:
  - "Async/sync wrapper pattern for Huey signal handlers"
  - "Optional notification field pattern (notification_email)"

# Metrics
duration: 4min
completed: 2026-01-19
---

# Phase 4 Plan 02: Email Notifications Summary

**Opt-in email notifications on job completion/failure using aiosmtplib with multipart HTML emails**

## Performance

- **Duration:** 4 min
- **Started:** 2026-01-19T21:30:49Z
- **Completed:** 2026-01-19T21:34:54Z
- **Tasks:** 3
- **Files modified:** 7

## Accomplishments

- Added notification_email field to JobInput and JobSubmitRequest (opt-in)
- Created notifications module with async send_job_notification and sync wrapper
- Integrated email sending into Huey signal handlers for completion and failure
- Multipart email with plain text and HTML versions, includes all download links

## Task Commits

Each task was committed atomically:

1. **Task 1: Add aiosmtplib dependency and notification email field** - `eef3f14` (feat)
2. **Task 2: Create notifications module** - `fcc16b3` (feat)
3. **Task 3: Update API and signal handlers** - `d2b7504` (feat)

## Files Created/Modified

- `src/qm_nmr_calc/notifications.py` - Email notification functions (async and sync)
- `pyproject.toml` - Added aiosmtplib and email-validator dependencies
- `src/qm_nmr_calc/models.py` - Added notification_email to JobInput
- `src/qm_nmr_calc/api/schemas.py` - Added notification_email with EmailStr validation
- `src/qm_nmr_calc/api/routers/jobs.py` - Pass notification_email in both endpoints
- `src/qm_nmr_calc/storage.py` - Accept notification_email in create_job_directory
- `src/qm_nmr_calc/queue.py` - Send emails in signal handlers

## Decisions Made

- **aiosmtplib over smtplib**: Async library for non-blocking email delivery
- **email-validator for EmailStr**: Required by Pydantic for email validation, returns 422 on invalid
- **Best-effort delivery**: Email failures are logged but never fail jobs
- **Environment variables for config**: SMTP_HOST, SMTP_PORT, SMTP_USER, SMTP_PASSWORD, NOTIFICATION_FROM, BASE_URL
- **Sync wrapper for Huey**: asyncio.run() wrapper since Huey signals are synchronous

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Missing email-validator dependency**
- **Found during:** Task 1 (verification)
- **Issue:** Pydantic EmailStr requires email-validator package, not installed
- **Fix:** Added email-validator to dependencies via `uv add email-validator`
- **Files modified:** pyproject.toml, uv.lock
- **Verification:** JobSubmitRequest with EmailStr imports successfully
- **Committed in:** eef3f14 (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Necessary dependency for EmailStr validation. No scope creep.

## Issues Encountered

None - plan executed smoothly after adding missing dependency.

## User Setup Required

**External services require manual configuration.** Environment variables needed for email delivery:

| Variable | Description | Example |
|----------|-------------|---------|
| SMTP_HOST | SMTP server hostname | smtp.gmail.com |
| SMTP_PORT | SMTP server port | 587 |
| SMTP_USER | SMTP username/API key | your-email@gmail.com |
| SMTP_PASSWORD | SMTP password/API key | app-specific-password |
| NOTIFICATION_FROM | Sender email address | noreply@qm-nmr-calc.example |
| BASE_URL | Public URL of the service | https://qm-nmr-calc.example.com |

## Next Phase Readiness

- Email notifications ready for Phase 5 (Monitoring) - can alert on job completion
- All results delivery features complete (downloads + notifications)
- Phase 4 complete after this plan

---
*Phase: 04-results-delivery*
*Completed: 2026-01-19*
