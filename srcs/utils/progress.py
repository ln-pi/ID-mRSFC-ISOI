"""Worker progress messages. Numerical routines stay silent outside a batch."""
import json
import os

PREFIX = '@@ISOI_PROGRESS@@'


def emit_progress(done=None, total=None, phase='', **extra):
    if os.environ.get('ISOI_MANAGED_WORKER') == '1':
        print(PREFIX + json.dumps(dict(done=done, total=total, phase=phase, **extra)), flush=True)


def progress_items(items, total=None, phase=''):
    """Count completed work, including units skipped or marked unavailable."""
    if total is None:
        total = len(items)
    emit_progress(0, total, phase)
    for done, item in enumerate(items, 1):
        yield item
        emit_progress(done, total, phase)


def console_input(prompt):
    if os.environ.get('ISOI_MANAGED_WORKER') != '1':
        return input(prompt)
    emit_progress(kind='prompt', text=prompt)
    try:
        return input()
    finally:
        emit_progress(kind='resume')
