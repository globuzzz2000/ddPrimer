#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Logging configuration module for the ddPrimer pipeline.
"""

import os
import sys
import logging
from tqdm import tqdm
from datetime import datetime

# Import package modules
from .config import Config


class ColoredFormatter(logging.Formatter):
    """Custom formatter that adds ANSI color codes for different log levels in debug mode."""
    
    # ANSI color codes
    COLORS = {
        'DEBUG': '\033[37m',     # White
        'INFO': '\033[36m',      # Cyan
        'WARNING': '\033[33m',   # Yellow
        'ERROR': '\033[31m',     # Red
        'CRITICAL': '\033[35m',  # Magenta
    }
    RESET = '\033[0m'  # Reset to normal
    
    def __init__(self, fmt=None, datefmt=None, use_colors=False):
        """
        Initialize the colored formatter.
        
        Args:
            fmt: Log format string
            datefmt: Date format string
            use_colors: Whether to use ANSI color codes
        """
        super().__init__(fmt, datefmt)
        self.use_colors = use_colors
    
    def format(self, record):
        """
        Format the log record with colors if enabled.
        
        Args:
            record: LogRecord instance
            
        Returns:
            str: Formatted log message
        """
        if self.use_colors and record.levelname in self.COLORS:
            # Get the original formatted message
            original_format = super().format(record)
            
            # Add color codes
            color_code = self.COLORS[record.levelname]
            return f"{color_code}{original_format}{self.RESET}"
        else:
            return super().format(record)


def setup_logging(debug: bool = False) -> str:
    """
    Configure logging for the application.
    
    Args:
        debug: Enable debug mode
        
    Returns:
        str: Path to the log file
    """
    # Use the Config.DEBUG_MODE if not explicitly set via command line
    debug_enabled = debug or Config.DEBUG_MODE
    
    log_level = logging.DEBUG if debug_enabled else logging.INFO
    
    # Different log formats based on debug mode
    if debug_enabled:
        log_format = '%(asctime)s - %(name)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(message)s'
        use_colors = True  # Enable colors in debug mode
    else:
        log_format = '%(message)s'  # Simpler format for regular use
        use_colors = False
    
    # Set up logging to file
    log_dir = os.path.join(os.path.expanduser("~"), ".ddPrimer", "logs")
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"ddPrimer_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    
    # Create file handler (always uses detailed format without colors)
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG if debug_enabled else logging.INFO)
    file_handler.setFormatter(logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(message)s'
    ))
    
    # Create console handler with appropriate format and colors
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    
    # Use ColoredFormatter for console output
    console_formatter = ColoredFormatter(log_format, use_colors=use_colors)
    console_handler.setFormatter(console_formatter)
    
    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG if debug_enabled else logging.INFO)
    
    # Clear existing handlers to avoid duplicates
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Add the handlers
    root_logger.addHandler(file_handler)
    root_logger.addHandler(console_handler)
    
    # Configure our specific logger
    logger = logging.getLogger("ddPrimer")
    
    if debug_enabled:
        logger.debug(f"Debug mode: {debug_enabled}")
        logger.debug(f"Log file: {log_file}")
        logger.debug(f"Python version: {sys.version}")
        logger.debug(f"Platform: {sys.platform}")
        logger.info("INFO messages will appear in bold in debug mode")
    else:
        logger.debug(f"Log file: {log_file}")
    
    # Set Config.DEBUG_MODE to ensure consistency
    Config.DEBUG_MODE = debug_enabled
    
    return log_file


class ProgressReporter:
    """Base class for progress reporting."""
    
    def start(self, total: int, description: str) -> None:
        """
        Start the progress reporter.
        
        Args:
            total: Total number of items to process
            description: Description of the current task
        """
        pass
    
    def update(self, amount: int = 1) -> None:
        """
        Update the progress by a given amount.
        
        Args:
            amount: Amount to increment the progress by
        """
        pass
    
    def finish(self) -> None:
        """Close the progress reporter."""
        pass


class TqdmProgressReporter(ProgressReporter):
    """Progress reporter implementation using tqdm."""
    
    def __init__(self):
        """Initialize the tqdm progress reporter."""
        self.pbar = None
    
    def start(self, total: int, description: str) -> None:
        """
        Start the progress bar.
        
        Args:
            total: Total number of items to process
            description: Description of the current task
        """
        self.pbar = tqdm(total=total, desc=description)
    
    def update(self, amount: int = 1) -> None:
        """
        Update the progress by a given amount.
        
        Args:
            amount: Amount to increment the progress by
        """
        if self.pbar:
            self.pbar.update(amount)
    
    def finish(self) -> None:
        """Close the progress bar."""
        if self.pbar:
            self.pbar.close()
            self.pbar = None


class LoggingProgressReporter(ProgressReporter):
    """Progress reporter implementation using logging."""
    
    def __init__(self, logger_name: str = "ddPrimer", log_interval: int = 10):
        """
        Initialize the logging progress reporter.
        
        Args:
            logger_name: Name of the logger to use
            log_interval: Interval (in percentage) for progress log messages
        """
        self.logger = logging.getLogger(logger_name)
        self.log_interval = log_interval
        self.total = 0
        self.current = 0
        self.description = ""
        self.last_percentage = 0
    
    def start(self, total: int, description: str) -> None:
        """
        Start the progress reporter.
        
        Args:
            total: Total number of items to process
            description: Description of the current task
        """
        self.total = total
        self.current = 0
        self.description = description
        self.last_percentage = 0
        self.logger.info(f"Starting {description} (0/{total}, 0%)")
    
    def update(self, amount: int = 1) -> None:
        """
        Update the progress by a given amount.
        
        Args:
            amount: Amount to increment the progress by
        """
        self.current += amount
        
        # Calculate current percentage
        if self.total > 0:
            percentage = int((self.current / self.total) * 100)
            
            # Log if we've passed a log interval threshold
            if percentage >= self.last_percentage + self.log_interval or self.current == self.total:
                self.logger.info(f"{self.description}: {self.current}/{self.total} ({percentage}%)")
                self.last_percentage = percentage
    
    def finish(self) -> None:
        """Finish the progress reporting."""
        if self.current < self.total:
            self.current = self.total
            self.logger.info(f"Completed {self.description}: {self.current}/{self.total} (100%)")


def get_progress_reporter() -> ProgressReporter:
    """
    Get an appropriate progress reporter based on configuration.
    
    Returns:
        ProgressReporter: A progress reporter instance
    """
    # Check if we should use a progress reporter at all
    if not Config.SHOW_PROGRESS:
        return ProgressReporter()  # Return the base class which does nothing
    
    # Try to use tqdm if available
    try:
        return TqdmProgressReporter()
    except ImportError:
        # Fall back to logging-based reporter
        return LoggingProgressReporter()