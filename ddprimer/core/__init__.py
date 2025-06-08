"""
Core processing modules for the ddPrimer pipeline.

This subpackage contains the main processing functionality:
- sequence_processor: Sequence filtering and processing
- primer_processor: Primer filtering and processing
- annotation_processor: GFF annotation processing
- primer3_processor: Interface with Primer3 for primer design
- blast_processor: BLAST specificity checking
- vienna_processor: Thermodynamic calculations using ViennaRNA
- vcf_processor: Variant calling and masking
"""

# Define explicit exports that are always available
__all__ = [
    'SequenceProcessor',
    'PrimerProcessor',
    'AnnotationProcessor',
    'SNPMaskingProcessor',
    'Primer3Processor',
    'BlastProcessor',
    'ViennaRNAProcessor'
]

# Core processors that are always available
from .sequence_processor import SequenceProcessor
from .primer_processor import PrimerProcessor
from .annotation_processor import AnnotationProcessor
from .SNP_masking_processor import SNPMaskingProcessor
from .primer3_processor import Primer3Processor
from .blast_processor import BlastProcessor
from .vienna_processor import ViennaRNAProcessor
