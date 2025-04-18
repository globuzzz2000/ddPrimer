#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline mode modules for ddPrimer.

This package contains the implementation of different operation modes:
- Standard mode: Uses FASTA, VCF, and GFF files
- Direct mode: Uses sequence data from CSV/Excel files
- MAF mode: Cross-species primer design using MAF alignments
"""

from . import common
from .alignment_mode import run as run_alignment_mode
from .direct_mode import run as run_direct_mode
from .standard_mode import run as run_standard_mode

__all__ = [
    'common',
    'standard_mode',
    'direct_mode',
    'alignment_mode',
    'run_standard_mode',
    'run_direct_mode',
    'run_alignment_mode'
]