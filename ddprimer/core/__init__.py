"""
Core processing modules for the ddPrimer pipeline.

This subpackage contains the main processing functionality:
- sequence_processor: Sequence filtering and processing
- primer_processor: Primer filtering and processing
- annotation_processor: GFF annotation processing
- primer3_processor: Interface with Primer3 for primer design
- blast_processor: BLAST specificity checking
- nupack_processor: Thermodynamic calculations (if NUPACK is installed)
- vienna_processor: Alternative thermodynamic calculations using ViennaRNA (if ViennaRNA is installed)
- thermo_processor: Flexible thermodynamic processor (can use either NUPACK or ViennaRNA)
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
    'ThermoProcessor'
]

# Core processors that are always available
from .sequence_processor import SequenceProcessor
from .primer_processor import PrimerProcessor
from .annotation_processor import AnnotationProcessor
from .SNP_masking_processor import SNPMaskingProcessor
from .primer3_processor import Primer3Processor
from .blast_processor import BlastProcessor
from .thermo_processor import ThermoProcessor  # This is always available

# Conditionally import NupackProcessor
try:
    import nupack
    from .nupack_processor import NupackProcessor
    __all__.append('NupackProcessor')
except ImportError:
    # NUPACK is not installed, skip the import
    pass

# Conditionally import ViennaRNAProcessor
try:
    import RNA
    from .vienna_processor import ViennaRNAProcessor
    __all__.append('ViennaRNAProcessor')
except ImportError:
    # ViennaRNA is not installed, skip the import
    pass