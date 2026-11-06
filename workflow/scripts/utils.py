# workflow/scripts/utils.py

"""
Shared logging configuration and OS-level stream redirection for workflow scripts.

This module provides core infrastructure utility functions and specialized handlers to enable
consistent logging, dynamic stream evaluation, secure file descriptor redirection (`os.dup2`),
and centralized Python warning interception across both standard workflow execution and
multiprocessing worker pools.
"""

import logging
import os
import sys
import warnings
from pathlib import Path
from typing import Any


class DynamicStreamHandler(logging.StreamHandler):
    """
    A specialized StreamHandler that dynamically evaluates `sys.stdout` upon log emission and flushing.

    Standard Python StreamHandler instances cache the `sys.stdout` file object at initialization time.
    When standard streams are later redirected at the file descriptor level (e.g., via `os.dup2` in
    `init_snakemake_log`), cached stream references become stale. This class dynamically updates its
    target stream to current `sys.stdout` before every emit and flush operation, ensuring logs are
    routed correctly to the active Snakemake log file.
    """

    def emit(self, record):
        self.stream = sys.stdout
        super().emit(record)

    def flush(self):
        self.stream = sys.stdout
        super().flush()


def init_snakemake_log(snakemake_obj: Any, base_dir: Path | None = None) -> None:
    """
    Redirects standard output (FD 1) and standard error (FD 2) at the OS level to a Snakemake log file.

    Performs low-level file descriptor duplication (`os.dup2`) to ensure that all output—including
    C-level subprocess streams, third-party C extensions, and Python standard streams—is securely captured
    within the target log file specified by Snakemake. Implements path traversal validation and strict
    file opening flags (`O_NOFOLLOW`, `0o600`) to prevent symlink exploits and illegal path escapes.

    Args:
        snakemake_obj (:obj:`Any`): The global `snakemake` object injected by the Snakemake workflow execution engine.
        base_dir (:obj:`pathlib.Path`, optional): Permitted base directory boundary for log file validation.
            Defaults to the current working directory (`Path.cwd()`).

    Raises:
        ValueError: If the target log file path resolves outside the permitted base directory boundaries.
        OSError: If opening or duplicating the low-level log file descriptor fails.
    """
    if not (hasattr(snakemake_obj, "log") and snakemake_obj.log):
        return

    effective_base = (base_dir or Path.cwd()).resolve()
    raw_log_path = Path(str(snakemake_obj.log[0]))
    log_path = raw_log_path.resolve()

    if not log_path.is_relative_to(effective_base):
        raise ValueError(f"Log path is invalid or outside the permitted boundaries: {log_path}")

    log_path.parent.mkdir(parents=True, exist_ok=True)

    sys.stdout.flush()
    sys.stderr.flush()

    open_flags = os.O_WRONLY | os.O_CREAT | os.O_TRUNC
    if hasattr(os, "O_NOFOLLOW"):
        open_flags |= os.O_NOFOLLOW

    try:
        log_fd = os.open(log_path, open_flags, 0o600)
    except OSError as err:
        raise OSError(f"Unable to open the log file securely: {err}") from err

    try:
        os.dup2(log_fd, 1)
        os.dup2(log_fd, 2)
    finally:
        os.close(log_fd)

    sys.stdout = open(1, "w", buffering=1, encoding="utf-8")
    sys.stderr = open(2, "w", buffering=1, encoding="utf-8")


def setup_logger(name: str = "workflow") -> logging.Logger:
    """
    Configures and returns a standardized, thread-safe Logger instance.

    Initializes a named `logging.Logger` configured with an `INFO` logging threshold, an explicit timestamp
    date format (`YYYY-MM-DD HH:MM:SS`), and a `DynamicStreamHandler` to guarantee log capture integrity
    following stream redirection. Prevents duplicate handler registration on re-entrant calls.

    Args:
        name (:obj:`str`, optional): The logger hierarchy name (default: "workflow").

    Returns:
        logger (:obj:`logging.Logger`): The fully configured logger instance.
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)

    if not logger.handlers:
        handler = DynamicStreamHandler()
        formatter = logging.Formatter(
            fmt="[%(asctime)s] %(levelname)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    return logger


def setup_warnings() -> None:
    """
    Configures centralized Python warning capture and formatting for workflow execution.

    Overrides `warnings.formatwarning` with a single-line linear format, enables global `UserWarning`
    capture, and routes all Python standard warnings into the dedicated `"py.warnings"` logger.
    Disables propagation to the root logger (`propagate = False`) to prevent duplicate log entries while
    ensuring warning events are formatted with timestamps and attached to `DynamicStreamHandler`.
    Acts as a worker initializer in `multiprocessing.Pool` instances to propagate warning configuration.
    """
    warnings.formatwarning = lambda msg, cat, fname, lineno, line=None: f"{fname}:{lineno}: {cat.__name__}: {msg}"
    warnings.simplefilter("always", UserWarning)
    logging.captureWarnings(True)

    warn_logger = logging.getLogger("py.warnings")
    warn_logger.setLevel(logging.WARNING)
    warn_logger.propagate = False

    if not warn_logger.handlers:
        handler = DynamicStreamHandler()
        formatter = logging.Formatter(
            fmt="[%(asctime)s] %(levelname)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        handler.setFormatter(formatter)
        warn_logger.addHandler(handler)
