# Phase 65: Documentation - Research

**Researched:** 2026-02-11
**Domain:** Technical documentation for scientific software
**Confidence:** HIGH

## Summary

Phase 65 is a documentation update phase that must bring SCALING-FACTORS.md and README.md into alignment with the newly expanded 13-solvent system completed in Phase 64. This is a pattern-repetition phase - Phase 58 (v2.8 milestone) performed an identical task when expanding from 3 to 7 solvents, providing a proven template.

The research reveals that documentation synchronization in scientific software requires strict value verification (JSON source → markdown tables), consistent formatting patterns, and cross-reference validation. The project already has 26 scaling factors (13 solvents × 2 nuclei) stored in `scaling_factors.json`, but documentation lags behind: SCALING-FACTORS.md has only 24 entries (missing vacuum), and README lists only 7 solvents.

This is a low-risk, straightforward documentation task with zero code changes. The Phase 58 plan provides an executable template that needs minor adaptation for the scale change (7→13 solvents, 14→26 factor sets).

**Primary recommendation:** Follow the Phase 58 pattern exactly - update SCALING-FACTORS.md main table to include all 26 factor sets with values extracted from scaling_factors.json, update README feature line and solvent table to list all 13 solvents, and verify all values match with automated script.

## Standard Stack

This is a pure documentation task requiring no libraries. The project uses standard markdown for all documentation.

### Core
| Tool | Version | Purpose | Why Standard |
|------|---------|---------|--------------|
| Markdown | CommonMark | Documentation format | Universal, readable, version-controllable |
| Python | 3.11+ | Verification scripts | Already project dependency |
| JSON | - | Source of truth for values | Project's data storage format |

### Supporting
| Tool | Version | Purpose | When to Use |
|------|---------|---------|-------------|
| grep | system | Value verification | Quick checks during editing |
| pytest | (project) | Test suite validation | Ensure no regressions |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Markdown tables | reStructuredText | RST more complex, no benefit here |
| Manual verification | No verification | Dangerous - errors would ship to users |
| Inline values | Reference to JSON | Less readable for users |

**Installation:**
No additional installation needed. All tools already in project environment.

## Architecture Patterns

### Recommended Documentation Structure

The project already has the correct structure:

```
qm-nmr-calc/
├── README.md                              # User-facing overview
├── docs/
│   ├── science.md                        # NMR methodology
│   ├── libraries.md                      # Integration details
│   └── ...                               # Other docs
├── data/benchmark/delta50/
│   └── SCALING-FACTORS.md                # Detailed factor documentation
└── src/qm_nmr_calc/data/
    └── scaling_factors.json              # Source of truth
```

### Pattern 1: Single Source of Truth
**What:** All numeric values for scaling factors exist in exactly one place: `scaling_factors.json`. Documentation files derive values from this source.

**When to use:** Always for scientific/calibration data that must be accurate.

**Example:**
```python
# Source: Phase 58 verification pattern
import json
factors = json.load(open('src/qm_nmr_calc/data/scaling_factors.json'))
slope = factors['B3LYP/6-311+G(2d,p)/1H/CHCl3']['slope']  # -0.9375324445748368
# Round for table: -0.9375 (4 decimal places)
```

### Pattern 2: Precision-Consistent Rounding
**What:** All values in SCALING-FACTORS.md table must use consistent precision:
- Slope: 4 decimal places
- Intercept: 2 decimal places
- R²: 4 decimal places
- MAE: 3 decimal places
- RMSD: 3 decimal places
- n: integer

**When to use:** When converting JSON float values to human-readable tables.

**Example:**
```python
# From scaling_factors.json
raw_slope = -0.9375324445748368
raw_intercept = 29.917572478615273
raw_r2 = 0.9952450408899208
raw_mae = 0.12419765454971371

# For SCALING-FACTORS.md table
table_slope = f"{raw_slope:.4f}"      # -0.9375
table_intercept = f"{raw_intercept:.2f}"  # 29.92
table_r2 = f"{raw_r2:.4f}"           # 0.9952
table_mae = f"{raw_mae:.3f}"         # 0.124
```

### Pattern 3: Alphabetical Solvent Ordering
**What:** Within each nucleus group (13C, 1H), solvents are listed alphabetically in SCALING-FACTORS.md table. This matches existing pattern and aids readability.

**Order:** Acetone, Acetonitrile, Benzene, CHCl3, DCM, DMF, DMSO, Methanol, Pyridine, THF, Toluene, Water (vacuum not included in original DELTA50, special case)

**Example:**
```markdown
| Nucleus | Solvent | Slope | Intercept | R^2 | MAE (ppm) | RMSD (ppm) | n |
|---------|---------|-------|-----------|-----|-----------|------------|---|
| 13C | Acetone | ... |
| 13C | Acetonitrile | ... |
| 13C | Benzene | ... |
...
| 1H | Acetone | ... |
| 1H | Acetonitrile | ... |
```

### Pattern 4: Automated Verification
**What:** After editing documentation, write a throwaway Python script that:
1. Reads scaling_factors.json
2. Parses SCALING-FACTORS.md table rows
3. Compares each value with appropriate rounding tolerance
4. Reports PASS/FAIL with specifics

**When to use:** Every time SCALING-FACTORS.md is edited with new values.

**Example:**
```python
# Source: Phase 58 Task 1 verification pattern
import json
import re

# Read source of truth
with open('src/qm_nmr_calc/data/scaling_factors.json') as f:
    factors = json.load(f)

# Read documentation
with open('data/benchmark/delta50/SCALING-FACTORS.md') as f:
    md_content = f.read()

# Extract table rows (regex or simple parsing)
rows = re.findall(r'\| (1H|13C) \| (\w+) \| ([-\d.]+)', md_content)

# Verify each row
for nucleus, solvent, slope_str in rows:
    key = f'B3LYP/6-311+G(2d,p)/{nucleus}/{solvent}'
    expected_slope = round(factors[key]['slope'], 4)
    actual_slope = round(float(slope_str), 4)
    assert expected_slope == actual_slope, f"Mismatch for {key}"

print("PASS: All values verified")
```

### Anti-Patterns to Avoid

- **Manual value transcription:** Never type numeric values from terminal output or paper notes into documentation. Always extract programmatically from scaling_factors.json to avoid transcription errors.

- **Inconsistent precision:** Don't mix "0.124" and "0.1244" in the same column. Pick a precision and stick to it for the entire table.

- **Updating docs without updating code:** This phase only updates docs because code was already updated in Phase 64. Never let docs and code diverge.

- **Skipping verification:** Don't commit updated SCALING-FACTORS.md without running an automated comparison against scaling_factors.json. Human proofreading is insufficient for 26 numeric entries.

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Markdown table generation | Manual string concatenation | Python tabulate or simple f-strings | Off-by-one errors in column alignment |
| JSON parsing | Custom text processing | `json.load()` | Handles edge cases, escaping |
| Floating point comparison | `actual == expected` | `abs(actual - expected) < tol` | Floating point rounding errors |
| Regex for table parsing | Complex multi-line regex | Line-by-line split + simple regex | Easier to debug when it fails |

**Key insight:** Even for "simple" documentation tasks, small automation prevents human error. A 50-line Python script can verify 156 numeric values (26 factor sets × 6 values each) faster and more accurately than manual review.

## Common Pitfalls

### Pitfall 1: Off-by-Two Error in Vacuum Placement
**What goes wrong:** SCALING-FACTORS.md currently has 24 factor sets. Adding vacuum requires inserting 2 rows (1H and 13C), not 1. Easy to forget one nucleus.

**Why it happens:** Thinking in "solvents" (13) rather than "factor sets" (26 = 13 × 2).

**How to avoid:** Requirements explicitly state "26 factor sets" and verification script should check `len(factor_sets) == 26`.

**Warning signs:** Grep for `vacuum` returns 1 match instead of 2, or test script reports 25 factor sets.

### Pitfall 2: Solvent Name Inconsistency
**What goes wrong:** Using different capitalization or spacing for solvent names between README, SCALING-FACTORS.md, and solvents.py. Examples: "chcl3" vs "CHCl3" vs "Chloroform" vs "chloroform".

**Why it happens:** Three different contexts need three different formats:
- Code keys (solvents.py): lowercase, e.g., "chcl3"
- Display names (README table): Title case with deuteration, e.g., "Chloroform (CDCl3)"
- SCALING-FACTORS.md: Mixed case matching NWChem convention, e.g., "CHCl3"

**How to avoid:** Use solvents.py SUPPORTED_SOLVENTS dict as source of truth for display names. Extract display name from parentheses for README. For SCALING-FACTORS.md, match the existing convention in the file (CHCl3, DMSO, not chcl3, dmso).

**Warning signs:** README lists "dcm" but SCALING-FACTORS.md has "DCM". Case mismatch.

### Pitfall 3: Confidence Interval Rounding Errors
**What goes wrong:** Confidence intervals in scaling_factors.json are arrays `[lower, upper]`. Forgetting to round both values to match table precision, or accidentally swapping them.

**Why it happens:** CI values are nested in JSON structure, easy to miss when extracting values.

**How to avoid:**
```python
ci_slope = factors[key]['ci_slope']
ci_lower = round(ci_slope[0], 4)  # 4 decimals for slope
ci_upper = round(ci_slope[1], 4)
table_ci = f"({ci_lower}, {ci_upper})"
```

**Warning signs:** CI values in table have wrong precision or reversed order (upper before lower).

### Pitfall 4: Cross-Reference Staleness
**What goes wrong:** docs/science.md and docs/libraries.md contain older statements like "Scaling factors are currently available for CHCl3, DMSO, and vacuum" which become false after expansion to 13 solvents.

**Why it happens:** Documentation files have redundant information that isn't centralized. Hard to find all mentions.

**How to avoid:** Grep entire docs/ directory for hardcoded solvent lists: `grep -r "CHCl3\|DMSO\|vacuum" docs/`. Update or remove stale claims. Consider replacing with "See README for full list" to avoid duplication.

**Warning signs:** Grep returns matches in docs/*.md that list old 3-solvent or 7-solvent sets.

### Pitfall 5: Precision Loss in JSON→Markdown Conversion
**What goes wrong:** Using wrong number of decimal places when rounding, leading to values that don't match the JSON source within tolerance.

**Why it happens:** Inconsistent rounding functions (round, floor, ceil) or wrong precision argument.

**How to avoid:** Standardize on Python's built-in `round()` with explicit decimal places:
- `round(value, 4)` for slope, r_squared
- `round(value, 2)` for intercept
- `round(value, 3)` for mae, rmsd
- `int(value)` for n_points

**Warning signs:** Verification script fails with "expected -0.9375, got -0.9376" due to rounding at wrong stage.

## Code Examples

Verified patterns for this documentation task:

### Extracting All Solvent Names from scaling_factors.json
```python
# Source: Project inspection
import json

with open('src/qm_nmr_calc/data/scaling_factors.json') as f:
    factors = json.load(f)

# Extract unique solvents
solvents = sorted(set(key.split('/')[-1] for key in factors.keys()))
print(f"Solvents: {solvents}")
# Output: ['Acetone', 'Acetonitrile', 'Benzene', 'CHCl3', 'DCM', 'DMF',
#          'DMSO', 'Methanol', 'Pyridine', 'THF', 'Toluene', 'Water', 'vacuum']
```

### Generating README Solvent Table from solvents.py
```python
# Source: Project pattern inspection
SUPPORTED_SOLVENTS = {
    "chcl3": "Chloroform (CDCl3)",
    "dmso": "Dimethylsulfoxide (DMSO-d6)",
    "vacuum": "Gas phase (no solvent)",
    # ... all 13 entries
}

# Generate table rows
for code, desc in SUPPORTED_SOLVENTS.items():
    # Use case column requires manual thought, not in dict
    print(f"| `{code}` | {desc} | [use case] |")
```

### Verification Script Skeleton
```python
# Source: Phase 58 Task 1 verify pattern
import json
import re

# Read JSON source of truth
with open('src/qm_nmr_calc/data/scaling_factors.json') as f:
    factors = json.load(f)

# Read markdown table
with open('data/benchmark/delta50/SCALING-FACTORS.md') as f:
    content = f.read()

# Extract table section (between "| Nucleus | Solvent" header and next "##")
table_start = content.find('| Nucleus | Solvent')
table_end = content.find('\n##', table_start)
table_section = content[table_start:table_end]

# Parse rows
rows = re.findall(
    r'\| (1H|13C) \| (\w+) \| ([-\d.]+) .* \| ([-\d.]+) .* \| ([-\d.]+) \| ([-\d.]+) \| ([-\d.]+) \| (\d+)',
    table_section
)

print(f"Found {len(rows)} factor sets in markdown")
assert len(rows) == 26, f"Expected 26 factor sets, found {len(rows)}"

# Verify each row
for nucleus, solvent, slope_str, intercept_str, r2_str, mae_str, rmsd_str, n_str in rows:
    key = f'B3LYP/6-311+G(2d,p)/{nucleus}/{solvent}'
    f = factors[key]

    # Compare with rounding tolerance
    assert abs(float(slope_str) - round(f['slope'], 4)) < 0.0001
    assert abs(float(intercept_str) - round(f['intercept'], 2)) < 0.01
    assert abs(float(r2_str) - round(f['r_squared'], 4)) < 0.0001
    assert abs(float(mae_str) - round(f['mae'], 3)) < 0.001
    assert abs(float(rmsd_str) - round(f['rmsd'], 3)) < 0.001
    assert int(n_str) == f['n_points']

print("PASS: All 26 factor sets verified")
```

### Counting Solvents in README Table
```bash
# Source: Standard grep pattern
# Count rows with backtick-quoted solvent codes
grep -c '^| `' README.md
# Should return 13 (one per solvent)
```

## State of the Art

This is a stable domain. Markdown documentation practices haven't changed significantly.

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Manual table editing | Script-generated + manual review | 2020s | Fewer transcription errors |
| Inline numeric values | Reference to data files | 2010s | Single source of truth |
| Word/PDF docs | Markdown in git | 2010s | Version control, diff-able |

**Deprecated/outdated:**
- **reStructuredText for simple docs**: Markdown has won for most projects unless using Sphinx for complex API docs
- **Wiki-based documentation**: Git-based markdown docs are now standard for open source projects

## Open Questions

Things that couldn't be fully resolved:

1. **Should vacuum appear first or last in solvent ordering?**
   - What we know: Current SCALING-FACTORS.md has vacuum missing. README places it third (after chcl3 and dmso). Alphabetically, "vacuum" would be last. But it's conceptually a baseline (no solvent).
   - What's unclear: Is there a semantic reason to place vacuum first or keep it in README's current position (third)?
   - Recommendation: For SCALING-FACTORS.md, insert vacuum at the END of each nucleus group (after Water) since it wasn't part of original DELTA50 solvents. For README, keep current position (third) to avoid gratuitous diff and because it's arbitrary.

2. **Are there cross-references in docs/ that need updating?**
   - What we know: Grep shows docs/science.md and docs/libraries.md mention old 3-solvent and 7-solvent sets. Phase 65 requirements say "documentation cross-references remain valid."
   - What's unclear: Are those statements within scope of Phase 65 or separate tasks?
   - Recommendation: Grep for hardcoded solvent lists and update any that are demonstrably wrong (e.g., "currently available for CHCl3, DMSO, vacuum" → "currently available for 13 solvents - see README"). Mark as verification step.

3. **What about the "Notes" section at end of SCALING-FACTORS.md?**
   - What we know: Phase 58 updated the Notes section to describe COSMO methodology. Current notes may reference 7 solvents or need updating.
   - What's unclear: Didn't read the full notes section (file is 1591 lines).
   - Recommendation: Read notes section, update any hardcoded counts or lists to match 13 solvents. Ensure methodology description remains accurate.

## Sources

### Primary (HIGH confidence)
- Project codebase inspection (scaling_factors.json, SCALING-FACTORS.md, README.md, solvents.py)
- Phase 58 documentation plan (58-01-PLAN.md) - proven pattern for identical task
- .planning/STATE.md - current system state (13 solvents, 26 factor sets)
- .planning/ROADMAP.md - Phase 65 requirements

### Secondary (MEDIUM confidence)
- [Markdown Best Practices](https://www.markdownguide.org/basic-syntax/) - Standard markdown formatting
- [Scientific Data Documentation](https://www.nature.com/sdata/submission-guidelines) - Table presentation standards
- [Technical Writing for Markdown](https://experienceleague.adobe.com/en/docs/contributor/contributor-guide/writing-essentials/markdown) - Adobe contributor guide

### Tertiary (LOW confidence)
- None. All findings verified with project codebase.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - No external libraries, pure markdown editing
- Architecture: HIGH - Phase 58 provides complete pattern, codebase inspection confirms structure
- Pitfalls: HIGH - Identified from Phase 58 task structure and common documentation errors

**Research date:** 2026-02-11
**Valid until:** 2026-03-11 (30 days - stable domain, markdown practices don't change rapidly)

**Research scope notes:**
- This is a pattern-repetition phase. Phase 58 (v2.8 milestone) performed the same task at smaller scale (7 solvents → 14 factor sets).
- Zero code changes required. Only markdown file updates.
- All source data (scaling_factors.json) already exists and validated in Phase 63-64.
- Main challenge is precision/accuracy (copying 156 numeric values correctly), solved by automated verification.
