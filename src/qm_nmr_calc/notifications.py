"""Email notification functions for job completion."""

import asyncio
import logging
import os
from email.message import EmailMessage
from typing import Optional

import aiosmtplib

logger = logging.getLogger(__name__)

# Configuration from environment (with sensible defaults for development)
SMTP_HOST = os.getenv("SMTP_HOST", "localhost")
SMTP_PORT = int(os.getenv("SMTP_PORT", "587"))
SMTP_USER = os.getenv("SMTP_USER", "")
SMTP_PASSWORD = os.getenv("SMTP_PASSWORD", "")
NOTIFICATION_FROM = os.getenv("NOTIFICATION_FROM", "noreply@qm-nmr-calc.example")
BASE_URL = os.getenv("BASE_URL", "http://localhost:8000")


async def send_job_notification(
    to_email: str,
    job_id: str,
    status: str,
    error_message: Optional[str] = None,
) -> bool:
    """Send job completion/failure notification email.

    Args:
        to_email: Recipient email address
        job_id: Job identifier
        status: Job status (complete, failed)
        error_message: Error details if failed

    Returns:
        True if sent successfully, False otherwise
    """
    msg = EmailMessage()
    msg["From"] = NOTIFICATION_FROM
    msg["To"] = to_email
    msg["Subject"] = f"NMR Calculation {status.title()}: {job_id}"

    results_url = f"{BASE_URL}/api/v1/jobs/{job_id}"

    # Plain text version
    if status == "complete":
        plain_body = f"""\
Your NMR calculation has completed successfully.

Job ID: {job_id}
Status: Complete

View your results: {results_url}

Available downloads:
- Chemical shifts (JSON): {results_url}/results
- Optimized geometry (XYZ): {results_url}/geometry
- Optimized geometry (SDF): {results_url}/geometry.sdf
- Raw NWChem output: {results_url}/output
"""
    else:
        plain_body = f"""\
Your NMR calculation has failed.

Job ID: {job_id}
Status: Failed
Error: {error_message or "Unknown error"}

View details: {results_url}
"""

    msg.set_content(plain_body)

    # HTML version
    if status == "complete":
        html_body = f"""\
<html>
<body style="font-family: Arial, sans-serif; line-height: 1.6; max-width: 600px;">
<h2 style="color: #2e7d32;">NMR Calculation Complete</h2>
<p>Your NMR calculation has completed successfully.</p>
<table style="border-collapse: collapse; margin: 15px 0;">
<tr><td style="padding: 5px 15px 5px 0;"><strong>Job ID:</strong></td><td>{job_id}</td></tr>
<tr><td style="padding: 5px 15px 5px 0;"><strong>Status:</strong></td><td style="color: #2e7d32;">Complete</td></tr>
</table>
<h3>Downloads</h3>
<ul>
<li><a href="{results_url}/results">Chemical Shifts (JSON)</a></li>
<li><a href="{results_url}/geometry">Optimized Geometry (XYZ)</a></li>
<li><a href="{results_url}/geometry.sdf">Optimized Geometry (SDF)</a></li>
<li><a href="{results_url}/output">Raw NWChem Output (ZIP)</a></li>
</ul>
<p style="margin-top: 20px;">
<a href="{results_url}" style="background-color: #2e7d32; color: white; padding: 10px 20px; text-decoration: none; border-radius: 5px; display: inline-block;">View Results</a>
</p>
</body>
</html>
"""
    else:
        html_body = f"""\
<html>
<body style="font-family: Arial, sans-serif; line-height: 1.6; max-width: 600px;">
<h2 style="color: #c62828;">NMR Calculation Failed</h2>
<p>Your NMR calculation encountered an error.</p>
<table style="border-collapse: collapse; margin: 15px 0;">
<tr><td style="padding: 5px 15px 5px 0;"><strong>Job ID:</strong></td><td>{job_id}</td></tr>
<tr><td style="padding: 5px 15px 5px 0;"><strong>Status:</strong></td><td style="color: #c62828;">Failed</td></tr>
<tr><td style="padding: 5px 15px 5px 0;"><strong>Error:</strong></td><td>{error_message or "Unknown error"}</td></tr>
</table>
<p><a href="{results_url}">View Details</a></p>
</body>
</html>
"""

    msg.add_alternative(html_body, subtype="html")

    try:
        # Only use authentication if credentials provided
        kwargs = {
            "hostname": SMTP_HOST,
            "port": SMTP_PORT,
        }
        if SMTP_USER and SMTP_PASSWORD:
            kwargs["username"] = SMTP_USER
            kwargs["password"] = SMTP_PASSWORD
            kwargs["start_tls"] = True

        await aiosmtplib.send(msg, **kwargs)
        logger.info(f"Sent notification email to {to_email} for job {job_id}")
        return True
    except Exception as e:
        # Log error but don't fail the job - email is best-effort
        logger.warning(f"Failed to send notification email to {to_email}: {e}")
        return False


def send_job_notification_sync(
    to_email: str,
    job_id: str,
    status: str,
    error_message: Optional[str] = None,
) -> bool:
    """Synchronous wrapper for use in Huey signal handlers.

    Huey signals run in sync context, so we need asyncio.run().
    """
    return asyncio.run(
        send_job_notification(to_email, job_id, status, error_message)
    )
