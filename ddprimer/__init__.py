#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ddPrimer: Droplet Digital PCR Primer Design Pipeline

A comprehensive and streamlined pipeline for designing primers and probes
specifically optimized for droplet digital PCR (ddPCR)
"""

import os
import sys

# Define version
__version__ = '1.0.0'

# Insert the current directory to the path for imports if running directly
if __name__ == '__main__':
    if os.path.dirname(__file__) not in sys.path:
        sys.path.insert(0, os.path.dirname(__file__))

# Make essential components available at package level
from .pipeline import run_pipeline, main
from .config import Config

__all__ = ['run_pipeline', 'main', 'Config', '__version__']