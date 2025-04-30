#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Logging configuration override for testing.

This module provides functionality to suppress excessive logging during testing.
Place this file in the tests directory.
"""

import logging
import sys


def configure_test_logging():
    """
    Configure logging settings to minimize output during tests.
    
    This function:
    1. Suppresses verbose logs from NUPACK
    2. Limits ddPrimer logs to warnings and above
    3. Redirects logs to a file if desired
    """
    # Configure root logger to minimal level
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.ERROR)
    
    # Configure ddPrimer logger to show only important messages during tests
    ddprimer_logger = logging.getLogger("ddPrimer")
    ddprimer_logger.setLevel(logging.WARNING)
    
    # Suppress NUPACK logging
    nupack_loggers = [
        logging.getLogger("nupack.rebind.render"),
        logging.getLogger("nupack"),
        logging.getLogger("nupack.design"),
        logging.getLogger("nupack.constants")
    ]
    for logger in nupack_loggers:
        logger.setLevel(logging.CRITICAL)
    
    # Set high log level for other noisy modules
    for module in ["matplotlib", "PIL", "tqdm"]:
        logging.getLogger(module).setLevel(logging.CRITICAL)


def redirect_logs_to_file(filename="test_log.txt"):
    """
    Redirect all logging output to a file instead of stdout.
    
    Args:
        filename (str): Path to the log file
    """
    root_logger = logging.getLogger()
    
    # Clear any existing handlers
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Add a file handler
    file_handler = logging.FileHandler(filename, mode="w")
    file_handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    root_logger.addHandler(file_handler)


# Automatically configure logging when this module is imported
configure_test_logging()