#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline mode modules for ddPrimer.

This package contains the implementation of different operation modes:
- Standard mode: Uses FASTA, VCF, and GFF files
- Direct mode: Uses sequence data from CSV/Excel files
- Alignment mode: Primer design using MAF alignments
"""

from . import common
from .direct_mode import run as run_direct_mode
from .standard_mode import run as run_standard_mode

__all__ = [
    'common',
    'standard_mode',
    'direct_mode',
    'run_standard_mode',
    'run_direct_mode'
]