#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility to reduce excessive logging in test scripts.
This code should be included at the top of each test script.
"""

import logging
import os

# Configure root logger to only show WARNING and above
logging.basicConfig(
    level=logging.WARNING,
    format='%(levelname)s: %(message)s'
)

# Get logger for the test script
logger = logging.getLogger("test")
# Set this logger to INFO level so our test messages still show
logger.setLevel(logging.INFO)

# Disable other verbose loggers
logging.getLogger("nupack").setLevel(logging.ERROR)
logging.getLogger("ddPrimer").setLevel(logging.WARNING)

# Also disable any handlers that might be already configured
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

# Configure a clean handler for our test logger
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
logger.addHandler(handler)

# Set environment variable to disable verbose NUPACK output
os.environ["NUPACK_DISABLE_LOGGING"] = "1"