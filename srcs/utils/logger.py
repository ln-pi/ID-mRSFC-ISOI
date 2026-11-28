import logging
import os
import sys

from datetime import datetime
from pathlib import Path
from logging.handlers import RotatingFileHandler
from typing import Optional, cast

from tqdm import tqdm


class CustomLogger(logging.Logger):
    def raise_(
            self,
            exc: Exception,
            msg: Optional[str] = None,
            level: int = logging.ERROR,
            exc_info: bool = True,
    ) -> None:
        if msg is None:
            msg = str(exc)
        self.log(level, msg, exc_info=exc_info)
        raise exc
    
    def exc(
            self,
            msg: str,
            *args, **kwargs
    ) -> None:
        self.error(msg, *args, exc_info=True, **kwargs)

logging.setLoggerClass(CustomLogger)


class LevelFormatter(logging.Formatter):
    def __init__(self, fmt_full, fmt_smpl, smpl_lvl=logging.INFO, datafmt=None):
        super().__init__(fmt_full, datafmt)
        self.fmt_full = fmt_full
        self.fmt_smpl = fmt_smpl
        self.smpl_lvl = smpl_lvl

        self.full_formatter = logging.Formatter(fmt_full, datafmt)
        self.smpl_formatter = logging.Formatter(fmt_smpl)
    
    def format(self, record):
        # Set level
        if record.levelno <= self.smpl_lvl:
            return self.smpl_formatter.format(record)

        return self.full_formatter.format(record)


class TqdmLoggingHandler(logging.Handler):
    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.write(msg, file=sys.stderr)
            self.flush()
        except Exception:
            self.handleError(record)


def setup_logger(
        log_file: Optional[str] = None,
        level: int = logging.INFO,
        to_console: bool = True
) -> CustomLogger:
    logger = cast(CustomLogger, logging.getLogger("project_logger"))

    full_fmt = "[%(levelname)s]%(asctime)s - %(module)s:%(lineno)d - %(message)s"
    smpl_fmt = "%(message)s"

    formatter = LevelFormatter(
        fmt_full=full_fmt,
        fmt_smpl=smpl_fmt,
        smpl_lvl=logging.INFO,
        datafmt="%Y-%m-%d %H:%M:%S"
    )

    # Imports and planning are console-only. The entry point supplies a file
    # once it has allocated a batch; library calls do not redirect that file.
    logger.setLevel(level)
    logger.propagate = False
    if to_console and not any(isinstance(h, TqdmLoggingHandler) for h in logger.handlers):
        console = TqdmLoggingHandler()
        console.setFormatter(formatter)
        logger.addHandler(console)
    if log_file is not None:
        requested_log_file = str(Path(log_file).resolve())
        if getattr(logger, "_log_file", None) != requested_log_file:
            close_log_file()
            Path(requested_log_file).parent.mkdir(parents=True, exist_ok=True)
            file_handler = RotatingFileHandler(
                requested_log_file, maxBytes=5*1024*1024, backupCount=3, encoding="utf-8"
            )
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
            logger._log_file = requested_log_file
            _write_run_header(requested_log_file)
    
    return logger


def close_log_file():
    logger = logging.getLogger("project_logger")
    for handler in list(logger.handlers):
        if isinstance(handler, RotatingFileHandler):
            logger.removeHandler(handler)
            handler.close()
    logger._log_file = None


def _write_run_header(log_file):
    header = (
        "\n\n"
        "=============== RUN START ===============\n"
        f"Timestamp: {datetime.now()}\n"
        f"Script: {os.path.basename(sys.argv[0])}\n"
        "-----------------------------------------\n"
    )
    with open(log_file, "a", encoding="utf-8") as f:
        f.write(header)
