# Phase 6: Web UI - Research

**Researched:** 2026-01-20
**Domain:** FastAPI HTML serving, CSS frameworks, JavaScript UI patterns
**Confidence:** HIGH

## Summary

This phase adds a browser-based interface wrapping the existing API. The established stack for FastAPI HTML apps is Jinja2 templates with StaticFiles mounting. For minimal CSS with professional styling, Pico CSS v2 is the optimal choice - it provides semantic HTML styling without requiring utility classes, includes a precompiled blue theme matching user requirements, and weighs under 10KB gzipped.

The web UI needs three pages (submit, status, results) with polling for status updates and a simple lightbox for image enlargement. All can be achieved with vanilla JavaScript and no external JS dependencies beyond what Pico CSS provides.

**Primary recommendation:** Use Jinja2 templates with Pico CSS blue theme via CDN, vanilla JavaScript polling with recursive setTimeout pattern, and inline modal implementation using HTML dialog element.

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Jinja2 | 3.1+ | Template engine | Built into FastAPI via Starlette, industry standard |
| Pico CSS | 2.1.1 | CSS framework | Semantic HTML, minimal classes, blue theme available |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| basicLightbox | 5.0.4 | Image modal | Optional - only if native dialog needs more features |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Pico CSS | Bootstrap 5 | Bootstrap has more components but 15x larger, requires classes |
| Pico CSS | Tailwind CDN | More flexible but requires utility classes, larger payload |
| Pico CSS | Plain CSS | Full control but significant development time for professional look |
| basicLightbox | Native dialog | Native dialog is sufficient for this use case, zero dependencies |

**Installation:**
```bash
# Jinja2 is already available via FastAPI/Starlette
pip install jinja2

# CSS and JS loaded via CDN - no Python packages needed
```

**CDN Links (add to templates):**
```html
<!-- Pico CSS Blue Theme -->
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/@picocss/pico@2/css/pico.blue.min.css">
```

## Architecture Patterns

### Recommended Project Structure
```
src/qm_nmr_calc/
├── api/
│   ├── app.py              # Add template/static mounting here
│   ├── routers/
│   │   ├── jobs.py         # Existing API routes
│   │   └── web.py          # NEW: Web UI routes (not /api/v1 prefix)
│   └── templates/           # Jinja2 templates
│       ├── base.html       # Layout template with header, CDN links
│       ├── submit.html     # Submission form
│       ├── status.html     # Status page with polling
│       └── results.html    # Results display
└── static/                  # Static files (CSS overrides, JS)
    ├── css/
    │   └── custom.css      # Minimal custom styles if needed
    └── js/
        └── polling.js      # Status polling logic
```

### Pattern 1: FastAPI Template Setup
**What:** Configure Jinja2Templates and StaticFiles in app.py
**When to use:** Always for HTML serving
**Example:**
```python
# Source: https://fastapi.tiangolo.com/advanced/templates/
from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from pathlib import Path

app = FastAPI(...)

# Get the directory where this file is located
BASE_DIR = Path(__file__).resolve().parent

# Mount static files
app.mount("/static", StaticFiles(directory=BASE_DIR / "static"), name="static")

# Configure templates
templates = Jinja2Templates(directory=BASE_DIR / "templates")
```

### Pattern 2: Template Response (FastAPI 0.108.0+)
**What:** Return HTML from route using TemplateResponse
**When to use:** All web UI routes
**Example:**
```python
# Source: https://fastapi.tiangolo.com/advanced/templates/
from fastapi import APIRouter, Request
from fastapi.responses import HTMLResponse

router = APIRouter(tags=["web"])

@router.get("/", response_class=HTMLResponse)
async def home(request: Request):
    return templates.TemplateResponse(
        request=request,
        name="submit.html",
        context={"solvents": get_supported_solvents()}
    )
```

### Pattern 3: Recursive setTimeout Polling
**What:** Poll API for status updates without request queue buildup
**When to use:** Status page auto-refresh
**Example:**
```javascript
// Source: https://dev.to/igadii/think-twice-before-using-setinterval-for-api-polling-it-might-not-be-ideal-3n3
const POLL_INTERVAL = 3000; // 3 seconds

async function checkStatus(jobId) {
    try {
        const response = await fetch(`/api/v1/jobs/${jobId}`);
        const data = await response.json();

        updateStatusDisplay(data);

        if (data.status === 'complete') {
            // Redirect to results page
            window.location.href = `/results/${jobId}`;
        } else if (data.status === 'failed') {
            // Stay on page, show error
            showError(data.error_message);
        } else {
            // Still running - schedule next poll
            setTimeout(() => checkStatus(jobId), POLL_INTERVAL);
        }
    } catch (error) {
        console.error('Polling error:', error);
        // Retry after longer delay on error
        setTimeout(() => checkStatus(jobId), POLL_INTERVAL * 2);
    }
}

// Start polling when page loads
document.addEventListener('DOMContentLoaded', () => {
    const jobId = document.body.dataset.jobId;
    checkStatus(jobId);
});
```

### Pattern 4: Native HTML Dialog for Image Modal
**What:** Use semantic dialog element for lightbox (Pico CSS styled)
**When to use:** Image enlargement on results page
**Example:**
```html
<!-- Source: https://picocss.com/docs/modal -->
<dialog id="image-modal">
    <article>
        <header>
            <button aria-label="Close" rel="prev" onclick="closeModal()"></button>
        </header>
        <img id="modal-image" src="" alt="Enlarged view">
    </article>
</dialog>

<script>
function openModal(imageSrc) {
    document.getElementById('modal-image').src = imageSrc;
    document.getElementById('image-modal').showModal();
    document.body.classList.add('modal-is-open');
}

function closeModal() {
    document.getElementById('image-modal').close();
    document.body.classList.remove('modal-is-open');
}
</script>
```

### Pattern 5: Pico CSS Form Styling
**What:** Semantic HTML forms with automatic styling
**When to use:** Submission form
**Example:**
```html
<!-- Source: https://picocss.com/docs/forms -->
<form action="/submit" method="post" enctype="multipart/form-data">
    <label>
        SMILES String
        <input type="text" name="smiles" placeholder="e.g., CCO for ethanol">
        <small>Enter a valid SMILES representation</small>
    </label>

    <label>
        Or Upload MOL/SDF File
        <input type="file" name="file" accept=".mol,.sdf">
    </label>

    <label>
        Solvent
        <select name="solvent" required>
            <option value="">Select a solvent...</option>
            <option value="chcl3">Chloroform (CDCl3)</option>
            <!-- populated from API -->
        </select>
    </label>

    <label>
        Preset
        <select name="preset">
            <option value="production">Production (accurate)</option>
            <option value="draft">Draft (fast)</option>
        </select>
    </label>

    <label>
        Email Notification (optional)
        <input type="email" name="notification_email" placeholder="your@email.com">
        <small>Leave blank if you don't want email notification</small>
    </label>

    <button type="submit">Calculate NMR</button>
</form>
```

### Anti-Patterns to Avoid
- **setInterval for polling:** Creates request queue buildup if server is slow; use recursive setTimeout instead
- **JavaScript frameworks:** Overkill for 3 static pages; adds complexity and bundle size
- **Custom CSS from scratch:** Time-consuming for professional look; Pico CSS handles this
- **Separate API calls from form submit:** Use standard form submission with redirect to status page

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Clean form styling | Custom CSS for inputs | Pico CSS semantic styling | Already handles validation states, file inputs, selects |
| Modal/lightbox | Custom overlay + positioning | HTML dialog element | Native, accessible, Pico CSS styles it |
| Responsive layout | Media queries from scratch | Pico CSS container + grid | Built-in responsive behavior |
| Button styling | Custom hover/focus states | Pico CSS button defaults | Handles all states including .secondary |
| Color theming | CSS variables from scratch | Pico CSS blue theme CDN | Pre-configured blue primary color |
| Loading indicator | Custom spinner | Pico CSS aria-busy on button | `<button aria-busy="true">Loading...</button>` |

**Key insight:** Pico CSS is designed for exactly this use case - semantic HTML with minimal classes. Fighting against it by adding custom classes defeats the purpose.

## Common Pitfalls

### Pitfall 1: Wrong TemplateResponse Syntax (FastAPI version mismatch)
**What goes wrong:** TypeError when using old TemplateResponse API
**Why it happens:** FastAPI 0.108.0 changed the API; old examples show `request` in context dict
**How to avoid:** Use named parameters: `templates.TemplateResponse(request=request, name="template.html", context={...})`
**Warning signs:** Error about missing positional argument or unexpected keyword argument

### Pitfall 2: Static Files 404
**What goes wrong:** CSS/JS files return 404
**Why it happens:** StaticFiles directory path is relative to wrong location
**How to avoid:** Use `Path(__file__).resolve().parent` to get absolute path to templates/static
**Warning signs:** Files exist but server returns 404; works in dev but not production

### Pitfall 3: Polling Never Stops
**What goes wrong:** Status page keeps polling even after completion
**Why it happens:** Missing terminal conditions in polling loop
**How to avoid:** Check for 'complete', 'failed', AND error responses before scheduling next poll
**Warning signs:** Network tab shows requests continuing after job finishes

### Pitfall 4: Form File Upload Missing
**What goes wrong:** File upload field is empty on server
**Why it happens:** Form missing `enctype="multipart/form-data"`
**How to avoid:** Always add enctype to forms with file inputs
**Warning signs:** `file` parameter is None on server despite user selecting file

### Pitfall 5: Modal Doesn't Close on Outside Click
**What goes wrong:** Clicking backdrop doesn't close modal
**Why it happens:** Dialog element needs explicit handling for backdrop clicks
**How to avoid:** Add click handler checking if click target is dialog itself (backdrop)
**Warning signs:** Only close button works, backdrop click does nothing

### Pitfall 6: SMILES/File Mutual Exclusion Not Handled
**What goes wrong:** User submits both SMILES and file, unclear which is used
**Why it happens:** Form allows both inputs but backend prioritizes one
**How to avoid:** Add JavaScript to disable SMILES input when file selected (and vice versa), or show clear instruction
**Warning signs:** User confusion about which input method is being used

## Code Examples

Verified patterns from official sources:

### Base Template with Pico CSS
```html
<!-- templates/base.html -->
<!DOCTYPE html>
<html lang="en" data-theme="light">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <meta name="color-scheme" content="light dark">
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/@picocss/pico@2/css/pico.blue.min.css">
    <link rel="stylesheet" href="{{ url_for('static', path='/css/custom.css') }}">
    <title>{% block title %}QM NMR Calculator{% endblock %}</title>
</head>
<body>
    <header class="container">
        <nav>
            <ul>
                <li><strong>QM NMR Calculator</strong></li>
            </ul>
            <ul>
                <li><a href="/">Submit Job</a></li>
            </ul>
        </nav>
    </header>

    <main class="container">
        {% block content %}{% endblock %}
    </main>

    <footer class="container">
        <small>Powered by ISiCLE and NWChem</small>
    </footer>

    {% block scripts %}{% endblock %}
</body>
</html>
```

### Status Page with Elapsed Time
```html
<!-- templates/status.html -->
{% extends "base.html" %}

{% block content %}
<article>
    <header>
        <h2>Job Status: {{ job.job_id }}</h2>
    </header>

    <div id="status-display">
        <p><strong>Status:</strong> <span id="job-status">{{ job.status }}</span></p>
        <p><strong>Elapsed:</strong> <span id="elapsed-time">--</span></p>
        {% if job.input.name %}
        <p><strong>Molecule:</strong> {{ job.input.name }}</p>
        {% endif %}
        <p><strong>SMILES:</strong> <code>{{ job.input.smiles }}</code></p>
    </div>

    <div id="error-display" hidden>
        <p role="alert" class="error"></p>
    </div>

    <footer>
        <small>Page refreshes automatically while job is running</small>
    </footer>
</article>
{% endblock %}

{% block scripts %}
<script>
    const JOB_ID = "{{ job.job_id }}";
    const CREATED_AT = new Date("{{ job.created_at }}");
    const POLL_INTERVAL = 3000;

    function updateElapsedTime() {
        const now = new Date();
        const elapsed = Math.floor((now - CREATED_AT) / 1000);
        const minutes = Math.floor(elapsed / 60);
        const seconds = elapsed % 60;
        document.getElementById('elapsed-time').textContent =
            `${minutes}m ${seconds}s`;
    }

    async function checkStatus() {
        try {
            const response = await fetch(`/api/v1/jobs/${JOB_ID}`);
            const data = await response.json();

            document.getElementById('job-status').textContent = data.status;
            updateElapsedTime();

            if (data.status === 'complete') {
                window.location.href = `/results/${JOB_ID}`;
            } else if (data.status === 'failed') {
                document.getElementById('error-display').hidden = false;
                document.querySelector('#error-display .error').textContent =
                    data.error_message || 'Calculation failed';
            } else {
                setTimeout(checkStatus, POLL_INTERVAL);
            }
        } catch (error) {
            console.error('Polling error:', error);
            setTimeout(checkStatus, POLL_INTERVAL * 2);
        }
    }

    // Start polling and elapsed time updates
    checkStatus();
    setInterval(updateElapsedTime, 1000);
</script>
{% endblock %}
```

### Results Page with Image Grid and Modal
```html
<!-- templates/results.html (key sections) -->
{% extends "base.html" %}

{% block content %}
<h2>Results: {{ job.input.name or job.job_id }}</h2>

<!-- Metadata Card -->
<article>
    <header>Calculation Details</header>
    <div class="grid">
        <div><strong>Preset:</strong> {{ job.preset }}</div>
        <div><strong>Solvent:</strong> {{ job.solvent }}</div>
        <div><strong>Functional:</strong> {{ results.functional }}</div>
        <div><strong>Basis Set:</strong> {{ results.basis_set }}</div>
    </div>
</article>

<!-- Image Grid -->
<article>
    <header>Visualizations</header>
    <div class="grid">
        <figure>
            <img src="/api/v1/jobs/{{ job.job_id }}/structure.png"
                 alt="Annotated structure"
                 onclick="openModal(this.src)">
            <figcaption>Annotated Structure</figcaption>
        </figure>
        <figure>
            <img src="/api/v1/jobs/{{ job.job_id }}/spectrum/1h.png"
                 alt="1H NMR spectrum"
                 onclick="openModal(this.src)">
            <figcaption>1H NMR Spectrum</figcaption>
        </figure>
        <figure>
            <img src="/api/v1/jobs/{{ job.job_id }}/spectrum/13c.png"
                 alt="13C NMR spectrum"
                 onclick="openModal(this.src)">
            <figcaption>13C NMR Spectrum</figcaption>
        </figure>
    </div>
</article>

<!-- Downloads -->
<article>
    <header>Downloads</header>
    <div class="grid">
        <a href="/api/v1/jobs/{{ job.job_id }}/geometry" role="button" class="secondary">
            Geometry (XYZ)
        </a>
        <a href="/api/v1/jobs/{{ job.job_id }}/geometry.sdf" role="button" class="secondary">
            Geometry (SDF)
        </a>
        <a href="/api/v1/jobs/{{ job.job_id }}/output" role="button" class="secondary">
            Raw Output (ZIP)
        </a>
    </div>
</article>

<!-- Image Modal -->
<dialog id="image-modal">
    <article>
        <header>
            <button aria-label="Close" rel="prev" onclick="closeModal()"></button>
        </header>
        <img id="modal-image" src="" alt="Enlarged view" style="max-width: 100%; max-height: 80vh;">
    </article>
</dialog>
{% endblock %}

{% block scripts %}
<script>
    const modal = document.getElementById('image-modal');

    function openModal(src) {
        document.getElementById('modal-image').src = src;
        modal.showModal();
        document.body.classList.add('modal-is-open');
    }

    function closeModal() {
        modal.close();
        document.body.classList.remove('modal-is-open');
    }

    // Close on backdrop click
    modal.addEventListener('click', (e) => {
        if (e.target === modal) closeModal();
    });

    // Close on Escape
    modal.addEventListener('cancel', closeModal);
</script>
{% endblock %}
```

## Codebase Integration Points

### Existing API Endpoints (reuse via fetch or links)
| Endpoint | Purpose | Web UI Usage |
|----------|---------|--------------|
| `POST /api/v1/jobs` | Submit SMILES | Form submission (JSON) |
| `POST /api/v1/jobs/upload` | Submit file | Form submission (multipart) |
| `GET /api/v1/jobs/{id}` | Get status | Status page polling |
| `GET /api/v1/jobs/{id}/results` | Get NMR results | Results page data |
| `GET /api/v1/jobs/{id}/spectrum/1h.png` | 1H spectrum image | Results page display |
| `GET /api/v1/jobs/{id}/spectrum/13c.png` | 13C spectrum image | Results page display |
| `GET /api/v1/jobs/{id}/structure.png` | Structure image | Results page display |
| `GET /api/v1/jobs/solvents` | List solvents | Form dropdown population |

### New Web Routes Needed
| Route | Method | Purpose |
|-------|--------|---------|
| `/` | GET | Submission form page |
| `/submit` | POST | Process form, redirect to status |
| `/status/{job_id}` | GET | Status page with polling |
| `/results/{job_id}` | GET | Results display page |

### Files to Modify
| File | Changes |
|------|---------|
| `src/qm_nmr_calc/api/app.py` | Add StaticFiles mount, Templates setup, include web router |

### Files to Create
| File | Purpose |
|------|---------|
| `src/qm_nmr_calc/api/routers/web.py` | Web UI routes |
| `src/qm_nmr_calc/api/templates/base.html` | Layout template |
| `src/qm_nmr_calc/api/templates/submit.html` | Submission form |
| `src/qm_nmr_calc/api/templates/status.html` | Status page |
| `src/qm_nmr_calc/api/templates/results.html` | Results page |
| `src/qm_nmr_calc/api/static/css/custom.css` | Minimal custom styles |

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Pass request in context dict | Named request parameter | FastAPI 0.108.0 / Starlette 0.29.0 | TemplateResponse API changed |
| setInterval polling | Recursive setTimeout | Ongoing best practice | Prevents request queue buildup |
| jQuery lightbox | Native dialog element | HTML5 dialog support | No JS library needed |
| Bootstrap/heavy frameworks | Semantic CSS (Pico) | 2020s trend | Smaller payload, cleaner HTML |

**Deprecated/outdated:**
- `templates.TemplateResponse("template.html", {"request": request})` - old syntax before FastAPI 0.108.0
- jQuery-based UI patterns - vanilla JS is sufficient for simple interactions

## Open Questions

Things that couldn't be fully resolved:

1. **Form submission handling approach**
   - What we know: Both SMILES and file endpoints exist
   - What's unclear: Should web form use JS fetch or traditional form POST?
   - Recommendation: Use traditional form POST to `/submit` endpoint that processes and redirects - simpler, no JS required for submission

2. **Mutual exclusion of SMILES vs file input**
   - What we know: User should submit one or the other
   - What's unclear: Should JavaScript disable one field when other is filled?
   - Recommendation: Add JS to clear/disable opposite field, with clear UI indication

## Sources

### Primary (HIGH confidence)
- [FastAPI Templates Documentation](https://fastapi.tiangolo.com/advanced/templates/) - Official template setup guide
- [Pico CSS Documentation](https://picocss.com/docs) - v2.1.1 current documentation
- [Pico CSS jsDelivr CDN](https://cdn.jsdelivr.net/npm/@picocss/pico@2/css/) - Available CSS files including blue theme

### Secondary (MEDIUM confidence)
- [DEV.to: API Polling with setTimeout](https://dev.to/igadii/think-twice-before-using-setinterval-for-api-polling-it-might-not-be-ideal-3n3) - Verified polling pattern
- [Real Python: FastAPI Templates](https://realpython.com/fastapi-jinja2-template/) - Comprehensive tutorial

### Tertiary (LOW confidence)
- None - all findings verified with official sources

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Context7/official docs for FastAPI, official docs for Pico CSS
- Architecture: HIGH - Official FastAPI patterns, verified polling approach
- Pitfalls: HIGH - Common issues documented in official sources and verified tutorials

**Research date:** 2026-01-20
**Valid until:** 2026-03-20 (60 days - stable technologies)
