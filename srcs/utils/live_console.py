"""Live stage display and one-way worker output capture into the native logger."""
import json
import logging
import os
from queue import Queue, Empty
import re
import shutil
import subprocess
import sys
from threading import Thread
import time
from types import SimpleNamespace

from srcs.utils.progress import PREFIX


class StageDisplay:
    def __init__(self, name, stream=None):
        self.name, self.stream = name, stream or sys.stdout
        self.tty = self.stream.isatty()
        self.start = time.monotonic()
        self.last = 0
        self.width = 0
        self.done = self.total = None
        self.phase = 'initializing'
        self.waiting = False
        self.tick(force=True)

    def clear(self):
        if self.tty:
            self.stream.write('\r' + ' ' * self.width + '\r')
            self.stream.flush()

    def update(self, event):
        kind = event.get('kind')
        if kind == 'prompt':
            self.clear()
            self.stream.write(event['text'])
            self.stream.flush()
            self.waiting = True
        elif kind == 'resume':
            self.waiting = False
            self.stream.write('\n')
        else:
            self.done, self.total = event.get('done'), event.get('total')
            self.phase = event.get('phase', '')
        self.tick()

    def tick(self, force=False):
        now = time.monotonic()
        if self.waiting or (not force and now - self.last < (.15 if self.tty else 10)):
            return
        self.last = now
        elapsed = int(now - self.start)
        if self.total is not None and self.total > 0 and self.done is not None:
            fraction = min(1., max(0., self.done / self.total))
            count = int(fraction * 20)
            bar = '#' * count + '-' * (20 - count)
            detail = f'{self.done}/{self.total}'
        else:
            pos = elapsed % 18
            bar = ' ' * pos + '===' + ' ' * (17 - pos)
            detail = ''
        line = f'{self.name} ... [{bar}] {detail} {self.phase} {elapsed}s'.rstrip()
        if self.tty:
            line = line[:max(20, shutil.get_terminal_size().columns - 1)]
            self.stream.write('\r' + line.ljust(self.width))
            self.width = len(line)
        else:
            self.stream.write(line + '\n')
        self.stream.flush()

    def finish(self, status):
        self.clear()
        self.stream.write(f'{self.name} ... {status}\n')
        self.stream.flush()


def execute_stage(command, *, cwd, env, logger, display):
    """Only the parent writes the log; drain output while refreshing the bar."""
    process = subprocess.Popen(command, cwd=cwd, env=env, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True, encoding='utf-8',
                               errors='replace', bufsize=1)
    events = Queue()

    def read_output():
        try:
            for line in process.stdout:
                events.put(line.rstrip('\r\n'))
        finally:
            events.put(None)

    reader = Thread(target=read_output, daemon=True)
    reader.start()
    warnings = 0
    try:
        while True:
            try:
                line = events.get(timeout=.15)
            except Empty:
                display.tick()
                continue
            if line is None:
                break
            if line.startswith(PREFIX):
                try:
                    event = json.loads(line[len(PREFIX):])
                    display.update(event)
                    if event.get('kind') == 'prompt':
                        logger.info('User input requested: %s', event['text'])
                    continue
                except (ValueError, KeyError, TypeError):
                    pass  # Preserve malformed worker output for diagnosis.
            if line:
                warnings += bool(re.search(r'^\[WARNING\]|\b\w*Warning:', line))
                record = re.match(r'^\[(DEBUG|INFO|WARNING|ERROR|CRITICAL)\]\s*(.*)', line)
                if record:
                    logger.log(getattr(logging, record[1]), record[2])
                else:
                    logger.info(line)
            display.tick()
        return SimpleNamespace(returncode=process.wait(), warning_count=warnings)
    except BaseException:
        # Ctrl+C normally reaches both console processes. Give workers time to
        # persist interruption markers before terminating a stuck process tree.
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            if os.name == 'nt':
                subprocess.run(['taskkill', '/PID', str(process.pid), '/T', '/F'],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            else:
                process.terminate()
            process.wait()
        raise
    finally:
        reader.join(timeout=1)
        process.stdout.close()
