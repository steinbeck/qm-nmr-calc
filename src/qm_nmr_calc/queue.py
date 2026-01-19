"""Huey task queue configuration and signal handlers."""
from datetime import datetime
import traceback
from pathlib import Path

from huey import SqliteHuey, signals

from .storage import update_job_status


# Huey instance with SQLite storage
# fsync=True ensures durability across crashes
huey = SqliteHuey(
    'qm-nmr-calc',
    filename='./data/huey.db',
    fsync=True
)


@huey.signal(signals.SIGNAL_EXECUTING)
def on_task_start(signal, task):
    """Update job status when task starts executing."""
    # Task args: (job_id,)
    if not task.args:
        return
    job_id = task.args[0]
    try:
        update_job_status(
            job_id,
            status='running',
            started_at=datetime.utcnow()
        )
    except Exception:
        # Don't let status update failure stop the task
        pass


@huey.signal(signals.SIGNAL_COMPLETE)
def on_task_complete(signal, task):
    """Update job status when task completes successfully."""
    if not task.args:
        return
    job_id = task.args[0]
    try:
        update_job_status(
            job_id,
            status='complete',
            completed_at=datetime.utcnow()
        )
    except Exception:
        pass


@huey.signal(signals.SIGNAL_ERROR)
def on_task_error(signal, task, exc=None):
    """Update job status when task fails with exception."""
    if not task.args:
        return
    job_id = task.args[0]
    try:
        update_job_status(
            job_id,
            status='failed',
            completed_at=datetime.utcnow(),
            error_message=str(exc) if exc else 'Unknown error',
            error_traceback=traceback.format_exc()
        )
    except Exception:
        pass


@huey.signal(signals.SIGNAL_INTERRUPTED)
def on_task_interrupted(signal, task):
    """Update job status when task is interrupted (graceful shutdown)."""
    if not task.args:
        return
    job_id = task.args[0]
    try:
        update_job_status(
            job_id,
            status='failed',
            completed_at=datetime.utcnow(),
            error_message='Task interrupted - process shutdown'
        )
    except Exception:
        pass
