#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File preparation module for ddPrimer pipeline.

Contains functionality for:
1. VCF file validation and preparation (bgzip, indexing, normalization)
2. FASTA file indexing and validation
3. INFO/AF field addition to VCF files when missing
4. Chromosome name mapping and harmonization across files
5. Interactive file preparation with user consent

This module ensures all input files are properly formatted and compatible
before proceeding with the ddPrimer pipeline, minimizing downstream errors
and improving processing reliability.
"""

import os
import subprocess
import tempfile
import gzip
import logging
import shutil
import re
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional, Set

# Import package modules
from ..config import Config, FileError, ExternalToolError, SequenceProcessingError
from ..config.logging_config import DebugLogLimiter

# Set up module logger
logger = logging.getLogger(__name__)


class FilePreparator:
    """
    Handles comprehensive file preparation and validation for ddPrimer pipeline.
    
    This class validates and prepares VCF, FASTA, and GFF files to ensure
    compatibility and proper formatting. It performs automatic fixes when
    possible and requests user permission for file modifications.
    """
    
    # Required external tools
    REQUIRED_TOOLS = [
        ('bcftools', 'bcftools --version'),
        ('bgzip', 'bgzip --version'),
        ('tabix', 'tabix --version'),
        ('samtools', 'samtools --version')
    ]
    
    def __init__(self):
        """Initialize file preparator with temporary directory and tool validation."""
        logger.debug("=== FILE PREPARATOR INITIALIZATION ===")
        
        try:
            self.temp_dir = os.path.join(Config.get_user_config_dir(), "temp")
            os.makedirs(self.temp_dir, exist_ok=True)
            self.prepared_files = {}
            
            self._validate_dependencies()
            logger.debug(f"Initialized FilePreparator with temp dir: {self.temp_dir}")
            
        except Exception as e:
            error_msg = f"Failed to initialize FilePreparator: {str(e)}"
            logger.error(error_msg)
            raise FileError(error_msg) from e
    
    #############################################################################
    #                           Main Workflow Methods
    #############################################################################
    
    @classmethod
    def prepare_pipeline_files_workflow(cls, vcf_file: Optional[str], fasta_file: str,
                                       gff_file: Optional[str] = None) -> Dict:
        """
        Prepare and validate all pipeline input files for workflow integration.
        This version is updated to gracefully handle optional VCF files.
        """
        logger.debug("=== WORKFLOW: FILE PREPARATION ===")
        preparator = None
        try:
            preparator = cls()
            preparator.set_reference_file(fasta_file)

            # Run the main preparation logic, which is now robust to None inputs
            result = preparator.prepare_files(vcf_file, fasta_file, gff_file)

            # If a VCF file was involved, generate the chromosome map for output formatting
            if vcf_file:
                # Use the potentially corrected VCF file path from the result
                processed_vcf = result.get('vcf_file', vcf_file)
                vcf_chroms = preparator._get_file_chromosomes(processed_vcf, 'vcf')
                fasta_seqs = preparator._get_file_chromosomes(fasta_file, 'fasta')
                chromosome_map = preparator._generate_standardized_chromosome_map(vcf_chroms, fasta_seqs)
                result['chromosome_map'] = chromosome_map
            else:
                result['chromosome_map'] = None  # No VCF, so no map needed

            logger.debug("=== END WORKFLOW: FILE PREPARATION ===")
            return result

        except Exception as e:
            logger.error(f"Error in file preparation workflow: {str(e)}")
            raise
        finally:
            if preparator:
                preparator.cleanup()
    
    def prepare_files(self, vcf_file: Optional[str], fasta_file: str, 
                     gff_file: Optional[str] = None) -> Dict:
        """
        Main entry point for file preparation workflow.
        
        Args:
            vcf_file: Optional path to VCF file
            fasta_file: Path to FASTA file  
            gff_file: Optional path to GFF file
            
        Returns:
            Dictionary with preparation results and final file paths
        """
        logger.debug("=== FILE PREPARATION WORKFLOW ===")
        logger.info("\nAnalyzing input files for compatibility...")
        
        try:
            # Validate input files exist (now handles optional vcf_file)
            self._validate_input_files(vcf_file, fasta_file, gff_file)
            self.set_reference_file(fasta_file)
            
            # Analyze files and determine what needs to be prepared
            analysis = self._analyze_all_files(vcf_file, fasta_file, gff_file)
            
            # Check if any preparation is needed
            if not analysis['needs_preparation']:
                logger.info("All files successfully validated")
                return self._create_success_result(vcf_file, fasta_file, gff_file, analysis, False)
            
            # Report issues and get user consent
            self._report_issues(analysis['issues'])
            if not self._get_user_consent(analysis['issues']):
                return self._create_failure_result('User declined file preparation', analysis['issues'])
            
            # Prepare files
            logger.info("\nCorrecting files...")
            prepared_files = self._prepare_all_files(vcf_file, fasta_file, gff_file, analysis)
            
            # Validate VCF preparation success if a VCF was provided
            if vcf_file and not self._validate_vcf_preparation(prepared_files, analysis):
                return self._create_failure_result('VCF file preparation failed', analysis['issues'])
            
            logger.info("File preparation completed successfully!\n")
            return self._create_success_result(
                prepared_files.get('vcf', vcf_file),
                prepared_files.get('fasta', fasta_file),
                prepared_files.get('gff', gff_file),
                analysis, True, prepared_files
            )
            
        except Exception as e:
            error_msg = f"Error in file preparation workflow: {str(e)}"
            logger.error(error_msg)
            return self._create_failure_result(error_msg, [])
    
    #############################################################################
    #                           File Analysis Methods
    #############################################################################
    
    def _analyze_all_files(self, vcf_file: Optional[str], fasta_file: str, 
                          gff_file: Optional[str] = None) -> Dict:
        """Analyze all files to determine what preparation is needed."""
        logger.debug("=== FILE ANALYSIS ===")
        
        issues = []
        
        # Analyze each file type, checking for None
        if vcf_file:
            issues.extend(self._analyze_single_file(vcf_file, 'vcf'))
        
        issues.extend(self._analyze_single_file(fasta_file, 'fasta'))

        if gff_file:
            issues.extend(self._analyze_single_file(gff_file, 'gff'))
        
        # Check chromosome compatibility only if VCF is provided
        if vcf_file:
            issues.extend(self._analyze_chromosome_compatibility(vcf_file, fasta_file, gff_file))
        
        logger.debug(f"Analysis complete: {len(issues)} issues found")
        return {
            'needs_preparation': len(issues) > 0,
            'issues': issues
        }
    
    def _analyze_single_file(self, file_path: str, file_type: str) -> List[Dict]:
        """Analyze a single file for issues based on its type."""
        issues = []
        
        if file_type == 'vcf':
            if not self._is_properly_bgzipped(file_path):
                issues.append(self._create_issue('vcf_compression', file_path, 
                                               'VCF file needs proper bgzip compression', 
                                               'Recompress with bgzip'))
            
            if not self._is_file_indexed(file_path, 'vcf'):
                issues.append(self._create_issue('vcf_index', file_path,
                                               'VCF file needs tabix indexing',
                                               'Create tabix index'))
            
            if not self._has_af_field(file_path):
                issues.append(self._create_issue('vcf_af_field', file_path,
                                               'VCF file lacks INFO/AF field',
                                               'Add INFO/AF field using bcftools +fill-tags'))
            
            if not self._is_vcf_normalized(file_path):
                issues.append(self._create_issue('vcf_normalization', file_path,
                                               'VCF file needs normalization',
                                               'Normalize variants with bcftools norm'))
        
        elif file_type == 'fasta':
            if not self._is_file_indexed(file_path, 'fasta'):
                issues.append(self._create_issue('fasta_index', file_path,
                                               'FASTA file needs samtools indexing',
                                               'Create samtools faidx index'))
        
        elif file_type == 'gff':
            if not file_path.endswith('.gz'):
                issues.append(self._create_issue('gff_index', file_path,
                                               'GFF file needs compression and tabix indexing',
                                               'Sort, compress and index GFF file'))
            elif not self._is_file_indexed(file_path, 'gff'):
                issues.append(self._create_issue('gff_index', file_path,
                                               'GFF file needs tabix indexing',
                                               'Create tabix index for compressed GFF'))
        
        return issues
    
    def _analyze_chromosome_compatibility(self, vcf_file: str, fasta_file: str,
                                        gff_file: Optional[str] = None) -> List[Dict]:
        """Analyze chromosome name compatibility between files."""
        issues = []
        
        try:
            # This method is only called if vcf_file is not None
            vcf_chroms = set(self._get_file_chromosomes(vcf_file, 'vcf').keys())
            fasta_chroms = set(self._get_file_chromosomes(fasta_file, 'fasta').keys())
            
            missing_in_fasta = vcf_chroms - fasta_chroms
            if missing_in_fasta:
                compat_analysis = self._check_chromosome_compatibility(vcf_file, fasta_file)
                issues.append({
                    'type': 'chromosome_mapping',
                    'files': [vcf_file, fasta_file],
                    'description': f'VCF chromosomes not found in FASTA: {", ".join(sorted(missing_in_fasta))}',
                    'action': 'Rename chromosomes in VCF to match FASTA',
                    'mapping': compat_analysis.get('suggested_mapping', {})
                })
            
            # Check GFF compatibility if provided
            if gff_file:
                gff_chroms = self._get_gff_chromosomes(gff_file)
                gff_only = set(gff_chroms) - fasta_chroms
                if gff_only:
                    issues.append(self._create_issue('gff_chromosome_mapping', 
                                                   [gff_file, fasta_file],
                                                   f'GFF contains chromosomes not in FASTA: {", ".join(sorted(gff_only))}',
                                                   'Filter GFF to match FASTA chromosomes'))
        
        except Exception as e:
            issues.append(self._create_issue('compatibility_check_failed', None,
                                           f'Could not verify chromosome compatibility: {str(e)}',
                                           'Manual verification recommended'))
        
        return issues
    
    #############################################################################
    #                           File Preparation Methods
    #############################################################################
    
    def _prepare_all_files(self, vcf_file: str, fasta_file: str,
                          gff_file: Optional[str], analysis: Dict) -> Dict:
        """Create corrected versions of files based on analysis."""
        prepared_files = {}
        
        # Group issues by file type for more efficient processing
        vcf_issues = [issue for issue in analysis['issues'] 
                     if issue.get('type', '').startswith('vcf') or 
                        issue.get('type') == 'chromosome_mapping']
        
        fasta_issues = [issue for issue in analysis['issues'] 
                       if issue.get('type', '').startswith('fasta')]
        
        gff_issues = [issue for issue in analysis['issues'] 
                     if issue.get('type', '').startswith('gff')] if gff_file else []
        
        # Prepare files as needed
        if vcf_issues:
            prepared_vcf = self._prepare_vcf_file(vcf_file, vcf_issues)
            if prepared_vcf:
                prepared_files['vcf'] = prepared_vcf
        
        if fasta_issues:
            prepared_fasta = self._prepare_fasta_file(fasta_file, fasta_issues)
            if prepared_fasta:
                prepared_files['fasta'] = prepared_fasta
        
        if gff_issues:
            prepared_gff = self._prepare_gff_file(gff_file, gff_issues)
            if prepared_gff:
                prepared_files['gff'] = prepared_gff
        
        return prepared_files
    
    def _prepare_vcf_file(self, vcf_file: str, issues: List[Dict]) -> Optional[str]:
        """Prepare corrected VCF file with streamlined processing."""
        logger.debug("=== VCF PREPARATION ===")
        
        try:
            output_path = self._create_output_path(vcf_file, 'prepared.vcf.gz')
            self._cleanup_existing_file(output_path)
            
            # Process VCF through pipeline stages
            current_file = vcf_file
            pipeline_stages = [
                ('compression', self._process_vcf_compression),
                ('af_field', self._process_vcf_af_field),
                ('chromosome_mapping', self._process_vcf_chromosome_mapping),
                ('normalization', self._process_vcf_normalization)
            ]
            
            for stage_name, processor in pipeline_stages:
                stage_issues = [issue for issue in issues if stage_name in issue.get('type', '')]
                if stage_issues:
                    new_file = processor(current_file, stage_issues)
                    if new_file and new_file != current_file:
                        self._cleanup_temp_file(current_file, vcf_file)
                        current_file = new_file
            
            # Copy to final location and create index
            shutil.copy2(current_file, output_path)
            self._cleanup_temp_file(current_file, vcf_file)
            self._index_file(output_path, 'vcf')
            
            # Final verification
            if not self._verify_file_integrity(output_path, 'vcf'):
                self._cleanup_existing_file(output_path)
                return None
            
            logger.info(f"Prepared VCF file: {output_path}")
            return output_path
            
        except Exception as e:
            logger.error(f"Error preparing VCF file: {str(e)}")
            return None
    
    #############################################################################
    #                           VCF Processing Pipeline
    #############################################################################
    
    def _process_vcf_compression(self, vcf_file: str, issues: List[Dict]) -> str:
        """Process VCF compression stage."""
        logger.debug("Processing VCF compression")
        temp_path = self._create_temp_file('.vcf.gz')
        
        if vcf_file.endswith('.gz') and self._is_properly_bgzipped(vcf_file):
            shutil.copy2(vcf_file, temp_path)
        else:
            self._run_command(['bgzip', '-c', vcf_file], output_file=temp_path)
        
        return temp_path
    
    def _process_vcf_af_field(self, vcf_file: str, issues: List[Dict]) -> str:
        """Process VCF AF field addition stage."""
        logger.debug("Adding INFO/AF field")
        temp_path = self._create_temp_file('.vcf.gz')
        
        cmd = ['bcftools', '+fill-tags', vcf_file, '-Oz', '-o', temp_path, '--', '-t', 'AF']
        self._run_command(cmd)
        
        return temp_path
    
    def _process_vcf_chromosome_mapping(self, vcf_file: str, issues: List[Dict]) -> str:
        """Process VCF chromosome mapping stage."""
        mapping_issue = next((issue for issue in issues if issue.get('type') == 'chromosome_mapping'), None)
        if not mapping_issue or not mapping_issue.get('mapping'):
            return vcf_file
        
        logger.debug("Applying chromosome mapping")
        temp_path = self._create_temp_file('.vcf.gz')
        rename_file = self._create_chromosome_rename_file(mapping_issue['mapping'])
        
        try:
            cmd = ['bcftools', 'annotate', vcf_file, '--rename-chrs', rename_file, '-Oz', '-o', temp_path]
            self._run_command(cmd)
            return temp_path
        finally:
            self._cleanup_temp_file(rename_file)
    
    def _process_vcf_normalization(self, vcf_file: str, issues: List[Dict]) -> str:
        """Process VCF normalization stage."""
        logger.debug("Normalizing VCF")
        temp_path = self._create_temp_file('.vcf.gz')
        
        reference_file = self._get_reference_file()
        if not os.path.exists(f"{reference_file}.fai"):
            self._index_file(reference_file, 'fasta')
        
        cmd = ['bcftools', 'norm', vcf_file, '-f', reference_file, '-m', '-any', '-N', '-Oz', '-o', temp_path]
        # Suppress bcftools norm statistics output by redirecting stderr to devnull
        self._run_command(cmd, suppress_stderr=True)
        
        return temp_path
    
    #############################################################################
    #                           File Processing Utilities
    #############################################################################
    
    def _get_file_chromosomes(self, file_path: str, file_type: str) -> Dict[str, int]:
        """Extract chromosome/sequence information from files."""
        if file_type == 'vcf':
            return self._get_vcf_chromosomes(file_path)
        elif file_type == 'fasta':
            return self._get_fasta_sequences(file_path)
        else:
            return {}
    
    def _get_vcf_chromosomes(self, vcf_file: str) -> Dict[str, int]:
        """Extract chromosome names and variant counts from VCF file."""
        try:
            result = self._run_command(['bcftools', 'query', '-f', '%CHROM\n', vcf_file], capture_output=True)
            chromosomes = result.stdout.strip().split('\n')
            chrom_counts = defaultdict(int)
            
            for chrom in chromosomes:
                if chrom and chrom.strip():
                    chrom_counts[chrom.strip()] += 1
            
            if DebugLogLimiter.should_log('vcf_chromosomes_extracted', interval=10, max_initial=2):
                logger.debug(f"Found {len(chrom_counts)} unique chromosomes in VCF")
            return dict(chrom_counts)
            
        except Exception as e:
            logger.error(f"Error reading VCF chromosomes: {str(e)}")
            raise ExternalToolError(f"Error reading VCF chromosomes: {str(e)}", tool_name="bcftools") from e
    
    def _get_fasta_sequences(self, fasta_file: str) -> Dict[str, int]:
        """Extract sequence names and lengths from FASTA file."""
        try:
            sequences = {}
            current_seq = None
            current_length = 0
            
            opener = gzip.open if fasta_file.endswith('.gz') else open
            mode = 'rt' if fasta_file.endswith('.gz') else 'r'
                
            with opener(fasta_file, mode) as f:
                for line in f:
                    line = line.strip()
                    
                    if line.startswith('>'):
                        if current_seq is not None:
                            sequences[current_seq] = current_length
                        
                        seq_id = line[1:].split()[0]
                        current_seq = seq_id
                        current_length = 0
                    elif line:
                        current_length += len(line)
                
                if current_seq is not None:
                    sequences[current_seq] = current_length
            
            if DebugLogLimiter.should_log('fasta_sequences_extracted', interval=10, max_initial=2):
                logger.debug(f"Found {len(sequences)} sequences in FASTA")
            return sequences
            
        except Exception as e:
            logger.error(f"Error reading FASTA sequences: {str(e)}")
            raise FileError(f"Error reading FASTA sequences: {str(e)}") from e
    
    def _filter_nuclear_chromosomes(self, fasta_seqs: Dict[str, int]) -> Dict[str, int]:
        """Filter FASTA sequences to keep only likely nuclear chromosomes."""
        nuclear_seqs = {}
        
        organellar_indicators = [
            'MT', 'MITO', 'MITOCHONDRIAL', 'MITOCHONDRION',
            'PT', 'PLASTID', 'CHLOROPLAST', 'CHLORO', 'PLASMID', 'PLAS'
        ]
        
        for seq_name, length in fasta_seqs.items():
            seq_upper = seq_name.upper()
            is_organellar = any(indicator in seq_upper for indicator in organellar_indicators)
            is_too_small = length < Config.MIN_CHROMOSOME_SIZE
            
            if not is_organellar and not is_too_small:
                nuclear_seqs[seq_name] = length
                if DebugLogLimiter.should_log('nuclear_sequence_included', interval=50, max_initial=3):
                    logger.debug(f"Including nuclear sequence: {seq_name} ({length:,} bp)")
            else:
                if is_organellar and DebugLogLimiter.should_log('organellar_excluded', interval=20, max_initial=2):
                    logger.debug(f"Excluding organellar sequence: {seq_name} ({length:,} bp)")
                elif is_too_small and DebugLogLimiter.should_log('small_sequence_excluded', interval=50, max_initial=3):
                    logger.debug(f"Excluding small sequence: {seq_name} ({length:,} bp)")
        
        logger.debug(f"Filtered to {len(nuclear_seqs)} nuclear sequences from {len(fasta_seqs)} total")
        return nuclear_seqs
    
    #############################################################################
    #                           Utility Methods
    #############################################################################
    
    def _run_command(self, cmd: List[str], capture_output: bool = False, 
                    output_file: Optional[str] = None, timeout: int = 300, 
                    suppress_stderr: bool = False) -> subprocess.CompletedProcess:
        """Run external command with proper error handling."""
        try:
            if output_file:
                with open(output_file, 'wb') as f:
                    stderr_target = subprocess.DEVNULL if suppress_stderr else subprocess.PIPE
                    result = subprocess.run(cmd, stdout=f, stderr=stderr_target, timeout=timeout)
            else:
                stderr_target = subprocess.DEVNULL if suppress_stderr else subprocess.PIPE
                stdout_target = subprocess.PIPE if capture_output else None
                result = subprocess.run(cmd, stdout=stdout_target, stderr=stderr_target, 
                                      text=True if capture_output else False, timeout=timeout)
            
            if result.returncode != 0:
                tool_name = cmd[0] if cmd else "unknown"
                # If stderr was suppressed, we can't show the error details
                error_details = result.stderr if hasattr(result, 'stderr') and result.stderr and not suppress_stderr else 'Command failed'
                error_msg = f"{tool_name} execution failed: {error_details}"
                raise ExternalToolError(error_msg, tool_name=tool_name)
            
            return result
            
        except subprocess.TimeoutExpired:
            tool_name = cmd[0] if cmd else "unknown"
            error_msg = f"{tool_name} execution timed out after {timeout} seconds"
            raise ExternalToolError(error_msg, tool_name=tool_name)
        except Exception as e:
            tool_name = cmd[0] if cmd else "unknown"
            error_msg = f"Error running {tool_name}: {str(e)}"
            raise ExternalToolError(error_msg, tool_name=tool_name) from e
    
    def _create_temp_file(self, suffix: str) -> str:
        """Create a temporary file path."""
        fd, temp_path = tempfile.mkstemp(suffix=suffix, prefix="ddprimer_", dir=self.temp_dir)
        os.close(fd)
        return temp_path
    
    def _create_output_path(self, input_file: str, suffix: str) -> str:
        """Create output path for prepared file."""
        base_name = Path(input_file).stem
        if base_name.endswith('.vcf'):
            base_name = base_name[:-4]
        
        output_dir = os.path.dirname(os.path.abspath(input_file))
        return os.path.join(output_dir, f"{base_name}_{suffix}")
    
    def _cleanup_temp_file(self, file_path: str, original_file: str = None):
        """Clean up temporary file if it's not the original."""
        if file_path and file_path != original_file and os.path.exists(file_path):
            try:
                os.unlink(file_path)
                for ext in ['.tbi', '.csi']:
                    idx_file = f"{file_path}{ext}"
                    if os.path.exists(idx_file):
                        os.unlink(idx_file)
            except OSError:
                pass
    
    def _cleanup_existing_file(self, file_path: str):
        """Remove existing file and its indices."""
        if os.path.exists(file_path):
            os.unlink(file_path)
            for ext in ['.tbi', '.csi']:
                idx_file = f"{file_path}{ext}"
                if os.path.exists(idx_file):
                    os.unlink(idx_file)
    
    def _create_issue(self, issue_type: str, file_path, description: str, action: str) -> Dict:
        """Create a standardized issue dictionary."""
        return {
            'type': issue_type,
            'file': file_path,
            'description': description,
            'action': action
        }
    
    def _create_success_result(self, vcf_file: str, fasta_file: str, gff_file: Optional[str],
                              analysis: Dict, changes_made: bool, prepared_files: Dict = None) -> Dict:
        """Create a success result dictionary."""
        return {
            'success': True,
            'vcf_file': vcf_file,
            'fasta_file': fasta_file,
            'gff_file': gff_file,
            'changes_made': changes_made,
            'issues_found': analysis['issues'],
            'prepared_files': prepared_files or {}
        }
    
    def _create_failure_result(self, reason: str, issues: List[Dict]) -> Dict:
        """Create a failure result dictionary."""
        return {
            'success': False,
            'reason': reason,
            'issues_found': issues
        }
    
    #############################################################################
    #                           Validation Methods
    #############################################################################
    
    def _validate_dependencies(self):
        """Validate that all required external tools are available."""
        missing_tools = []
        
        for tool_name, version_cmd in self.REQUIRED_TOOLS:
            try:
                result = subprocess.run(version_cmd.split(), capture_output=True, text=True, timeout=10)
                if result.returncode != 0:
                    missing_tools.append(tool_name)
            except (FileNotFoundError, subprocess.TimeoutExpired):
                missing_tools.append(tool_name)
        
        if missing_tools:
            error_msg = (
                f"Required tools not found: {', '.join(missing_tools)}\n"
                "Please install missing tools:\n"
                "  Ubuntu/Debian: sudo apt-get install bcftools bgzip tabix samtools\n"
                "  macOS: brew install bcftools htslib samtools\n"
                "  Conda: conda install -c bioconda bcftools htslib samtools"
            )
            raise ExternalToolError(error_msg, tool_name=missing_tools[0])
    
    def _validate_input_files(self, vcf_file: Optional[str], fasta_file: str, gff_file: Optional[str] = None):
        """Validate that input files exist and are accessible."""
        # Start with the always-required FASTA file
        files_to_check = [(fasta_file, "FASTA")]
        
        # Add optional files only if they are provided
        if vcf_file:
            files_to_check.append((vcf_file, "VCF"))
        if gff_file:
            files_to_check.append((gff_file, "GFF"))
        
        for file_path, file_type in files_to_check:
            if not os.path.exists(file_path):
                raise FileError(f"{file_type} file not found: {file_path}")
            if not os.access(file_path, os.R_OK):
                raise FileError(f"{file_type} file not readable: {file_path}")
    
    def _validate_vcf_preparation(self, prepared_files: Dict, analysis: Dict) -> bool:
        """Validate that VCF preparation actually succeeded."""
        vcf_issues = [issue for issue in analysis['issues'] 
                     if issue.get('type', '').startswith('vcf') or 
                        issue.get('type') == 'chromosome_mapping']
        
        if vcf_issues:
            vcf_path = prepared_files.get('vcf')
            if not vcf_path or not os.path.exists(vcf_path):
                logger.error("VCF file preparation failed - file not created")
                return False
        
        return True
    
    def _verify_file_integrity(self, file_path: str, file_type: str) -> bool:
        """Verify that a prepared file can be read properly."""
        try:
            if file_type == 'vcf':
                test_cmd = ['bcftools', 'view', '-h', file_path]
                result = subprocess.run(test_cmd, capture_output=True, text=True, timeout=30)
                return result.returncode == 0
            elif file_type == 'fasta':
                return os.path.exists(f"{file_path}.fai") or os.path.getsize(file_path) > 0
            return True
        except Exception:
            return False
    
    #############################################################################
    #                           File Status Check Methods
    #############################################################################
    
    def _is_properly_bgzipped(self, vcf_file: str) -> bool:
        """Check if VCF file is properly bgzip compressed."""
        try:
            result = subprocess.run(['bcftools', 'view', '-h', vcf_file], 
                                  capture_output=True, text=True, timeout=10)
            return result.returncode == 0
        except Exception:
            return False
    
    def _is_file_indexed(self, file_path: str, file_type: str) -> bool:
        """Check if file has appropriate index."""
        if file_type == 'vcf':
            return any(os.path.exists(f"{file_path}{ext}") for ext in ['.tbi', '.csi'])
        elif file_type == 'fasta':
            return os.path.exists(f"{file_path}.fai")
        elif file_type == 'gff':
            # For GFF files, we need to check if the file is compressed AND indexed
            # If it's not compressed, it needs preparation regardless of index status
            if not file_path.endswith('.gz'):
                return False  # Uncompressed GFF always needs preparation
            # If compressed, check for index
            return any(os.path.exists(f"{file_path}{ext}") for ext in ['.tbi', '.csi'])
        return False
    
    def _has_af_field(self, vcf_file: str) -> bool:
        """Check if VCF file has INFO/AF field."""
        try:
            result = subprocess.run(['bcftools', 'view', '-h', vcf_file], 
                                  capture_output=True, text=True, timeout=10)
            return result.returncode == 0 and '##INFO=<ID=AF,' in result.stdout
        except Exception:
            return False
    
    def _is_vcf_normalized(self, vcf_file: str) -> bool:
        """Check if VCF file appears to be normalized."""
        try:
            # Sample first 10000 variants for normalization check
            view_proc = subprocess.Popen(['bcftools', 'view', vcf_file], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            head_proc = subprocess.Popen(['head', '-n', '10000'], stdin=view_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            view_proc.stdout.close()
            
            norm_proc = subprocess.Popen(
                ['bcftools', 'norm', '--check-ref', 'w', '-f', self._get_reference_file(), '-'],
                stdin=head_proc.stdout, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True
            )
            head_proc.stdout.close()
            
            try:
                stdout, stderr = norm_proc.communicate(timeout=30)
                view_proc.wait()
                head_proc.wait()
                
                if norm_proc.returncode == 0:
                    # Check for normalization warnings
                    normalization_warnings = ['not left-aligned', 'not normalized', 'multiallelic']
                    has_warnings = any(warning in stderr.lower() for warning in normalization_warnings)
                    return not has_warnings
                else:
                    # If reference sequence errors, assume normalized
                    if "sequence not found" in stderr.lower():
                        return True
                    return False
                    
            except subprocess.TimeoutExpired:
                # Kill processes and assume normalized
                for proc in [norm_proc, head_proc, view_proc]:
                    try:
                        proc.kill()
                    except:
                        pass
                return True
                
        except Exception as e:
            if DebugLogLimiter.should_log('normalization_check_error', interval=10, max_initial=1):
                logger.debug(f"Error checking VCF normalization: {str(e)}")
            return False
    
    #############################################################################
    #                           File Processing Methods
    #############################################################################
    
    def _prepare_fasta_file(self, fasta_file: str, issues: List[Dict]) -> Optional[str]:
        """Prepare FASTA file (mainly indexing)."""
        try:
            index_issues = [issue for issue in issues if issue['type'] == 'fasta_index']
            if index_issues:
                self._index_file(fasta_file, 'fasta')
            return fasta_file
        except Exception as e:
            logger.error(f"Error preparing FASTA file: {str(e)}")
            return None
    
    def _prepare_gff_file(self, gff_file: str, issues: List[Dict]) -> Optional[str]:
        """Prepare GFF file (sorting and indexing)."""
        try:
            index_issues = [issue for issue in issues if issue['type'] == 'gff_index']
            if not index_issues:
                return gff_file
            
            # Check if file is already compressed
            if gff_file.endswith('.gz'):
                # File is compressed but needs indexing
                self._index_file(gff_file, 'gff')
                logger.info(f"Created index for compressed GFF file: {gff_file}")
                return gff_file
            else:
                # File needs compression and indexing
                output_path = self._create_output_path(gff_file, 'prepared.gff.gz')
                
                if os.path.exists(output_path):
                    user_input = input(f"\nPrepared GFF file exists: {output_path}\nOverwrite? [y/n]: ").strip().lower()
                    if user_input not in ['y', 'yes']:
                        # Make sure existing file is indexed
                        if not self._is_file_indexed(output_path, 'gff'):
                            self._index_file(output_path, 'gff')
                        return output_path
                    else:
                        self._cleanup_existing_file(output_path)
                
                self._sort_and_compress_gff(gff_file, output_path)
                self._index_file(output_path, 'gff')
                logger.info(f"Prepared GFF file: {output_path}")
                return output_path
            
        except Exception as e:
            logger.error(f"Error preparing GFF file: {str(e)}")
            return None
    
    def _index_file(self, file_path: str, file_type: str):
        """Create appropriate index for file type."""
        if file_type == 'vcf':
            self._run_command(['tabix', '-p', 'vcf', file_path])
        elif file_type == 'fasta':
            self._run_command(['samtools', 'faidx', file_path])
        elif file_type == 'gff':
            self._run_command(['tabix', '-p', 'gff', file_path])
    
    def _sort_and_compress_gff(self, gff_file: str, output_path: str):
        """Sort and compress GFF file for tabix indexing."""
        temp_sorted_path = self._create_temp_file('.gff')
        
        try:
            # Sort the GFF file
            sort_cmd = ['sort', '-k1,1', '-k4,4n', '-T', self.temp_dir, '-o', temp_sorted_path, gff_file]
            self._run_command(sort_cmd)
            
            # Compress with bgzip
            compress_cmd = ['bgzip', '-c', temp_sorted_path]
            self._run_command(compress_cmd, output_file=output_path)
            
        finally:
            self._cleanup_temp_file(temp_sorted_path)
    
    def _create_chromosome_rename_file(self, mapping: Dict[str, str]) -> str:
        """Create temporary file for chromosome renaming."""
        rename_file = self._create_temp_file('.txt')
        
        with open(rename_file, 'w') as f:
            for old_name, new_name in mapping.items():
                f.write(f"{old_name}\t{new_name}\n")
        
        return rename_file
    
    #############################################################################
    #                           Chromosome Analysis Methods
    #############################################################################
    
    def _check_chromosome_compatibility(self, vcf_file: str, fasta_file: str) -> Dict:
        """Check chromosome compatibility and suggest mapping."""
        try:
            vcf_chroms = self._get_file_chromosomes(vcf_file, 'vcf')
            fasta_seqs = self._get_file_chromosomes(fasta_file, 'fasta')
            
            vcf_set = set(vcf_chroms.keys())
            fasta_set = set(fasta_seqs.keys())
            
            exact_matches = vcf_set & fasta_set
            vcf_only = vcf_set - fasta_set
            fasta_only = fasta_set - vcf_set
            
            analysis = {
                'vcf_chromosomes': vcf_chroms,
                'fasta_sequences': fasta_seqs,
                'exact_matches': exact_matches,
                'vcf_only': vcf_only,
                'fasta_only': fasta_only,
                'compatible': len(exact_matches) > 0,
                'needs_mapping': len(vcf_only) > 0 or len(fasta_only) > 0
            }
            
            if analysis['needs_mapping']:
                analysis['suggested_mapping'] = self._suggest_chromosome_mapping(vcf_chroms, fasta_seqs)
            
            return analysis
            
        except Exception as e:
            logger.error(f"Error in compatibility check: {str(e)}")
            raise SequenceProcessingError(f"Error in compatibility check: {str(e)}") from e
    
    def _suggest_chromosome_mapping(self, vcf_chroms: Dict[str, int], 
                                   fasta_seqs: Dict[str, int]) -> Dict[str, str]:
        """Suggest intelligent chromosome mapping."""
        # Sort chromosomes by numeric component
        vcf_sorted = sorted(vcf_chroms.keys(), 
                           key=lambda x: (self._extract_numeric_component(x), x))
        
        nuclear_fasta = self._filter_nuclear_chromosomes(fasta_seqs)
        fasta_sorted = sorted(nuclear_fasta.keys(), 
                             key=lambda x: (self._extract_numeric_component(x), x))
        
        if len(vcf_sorted) <= len(fasta_sorted):
            main_fasta = fasta_sorted[:len(vcf_sorted)]
            mapping = dict(zip(vcf_sorted, main_fasta))
            return mapping
        
        return {}
    
    def _extract_numeric_component(self, name: str) -> float:
        """Extract numeric component from chromosome name for sorting."""
        if not name or not isinstance(name, str):
            return float('inf')
            
        name_upper = name.upper()
        
        # Simple numbers
        if name.isdigit():
            return int(name)
        
        # Chr1, Chr2, etc.
        if name_upper.startswith('CHR') and len(name) > 3:
            chr_part = name[3:]
            if chr_part.isdigit():
                return int(chr_part)
        
        # Extract numbers from accession-like names
        numbers = re.findall(r'\d+', name)
        if numbers:
            return int(numbers[-1])
        
        # Special chromosomes
        special_chroms = {
            'X': 100, 'Y': 101, 'MT': 102, 'MITO': 102, 'MITOCHONDRIAL': 102,
            'CHLOROPLAST': 103, 'PLASTID': 103, 'PT': 103, 'CP': 103
        }
        
        for special, value in special_chroms.items():
            if special in name_upper:
                return value
        
        return float('inf')
    
    def _generate_standardized_chromosome_map(self, vcf_chroms: Dict[str, int],
                                             fasta_seqs: Dict[str, int]) -> Dict[str, str]:
        """Generate mapping from original chromosome names to standardized names."""
        chrom_names = sorted(vcf_chroms.keys(),
                             key=lambda x: (self._extract_numeric_component(x), x))
        
        standardized_map = {}
        for i, name in enumerate(chrom_names):
            numeric_comp = self._extract_numeric_component(name)
            
            if numeric_comp == 100:
                standardized_name = 'X'
            elif numeric_comp == 101:
                standardized_name = 'Y'
            elif numeric_comp == 102:
                standardized_name = 'MT'
            elif numeric_comp < float('inf'):
                if name.isdigit():
                    standardized_name = name
                elif name.upper().startswith('CHR') and name[3:].isdigit():
                    standardized_name = name[3:]
                else:
                    standardized_name = str(i + 1)
            else:
                standardized_name = name
            
            standardized_map[name] = standardized_name
        
        return standardized_map
    
    def _get_gff_chromosomes(self, gff_file: str) -> Set[str]:
        """Extract chromosome names from GFF file."""
        chromosomes = set()
        
        try:
            opener = gzip.open if gff_file.endswith('.gz') else open
            mode = 'rt' if gff_file.endswith('.gz') else 'r'
            
            with opener(gff_file, mode) as f:
                for line_num, line in enumerate(f, 1):
                    if line_num > 1000:  # Sample first 1000 lines
                        break
                    
                    line = line.strip()
                    if line and not line.startswith('#'):
                        fields = line.split('\t')
                        if len(fields) >= 1:
                            chromosomes.add(fields[0])
        
        except Exception as e:
            if DebugLogLimiter.should_log('gff_chromosome_read_error', interval=10, max_initial=1):
                logger.debug(f"Error reading GFF chromosomes: {str(e)}")
        
        return chromosomes
    
    #############################################################################
    #                           User Interaction Methods
    #############################################################################
    
    def _report_issues(self, issues: List[Dict]):
        """Report found issues to user."""
        logger.info(f"Found {len(issues)} issue(s) that need to be addressed:\n")
        
        for i, issue in enumerate(issues, 1):
            logger.info(f"{i}. {issue['description']}")
    
    def _get_user_consent(self, issues: List[Dict]) -> bool:
        """Ask user for consent to create corrected files."""
        logger.info("\nTo proceed with the pipeline, corrected copies of your files need to be created.")
        
        while True:
            try:
                response = input("\n>>> Create corrected files? <<<\n[y/n]: ").strip().lower()
                if response in ['y', 'yes']:
                    return True
                elif response in ['n', 'no']:
                    return False
                else:
                    print("Please enter 'y' for yes or 'n' for no.")
            except (EOFError, KeyboardInterrupt):
                return False
    
    #############################################################################
    #                           Reference File Management
    #############################################################################
    
    def set_reference_file(self, reference_file: str):
        """Set reference FASTA file for operations that require it."""
        if not os.path.exists(reference_file):
            raise FileError(f"Reference FASTA file not found: {reference_file}")
        
        self._reference_file = reference_file
        if DebugLogLimiter.should_log('reference_file_set', interval=1, max_initial=1):
            logger.debug(f"Set reference file: {reference_file}")
    
    def _get_reference_file(self) -> str:
        """Get reference FASTA file for normalization."""
        if hasattr(self, '_reference_file') and self._reference_file:
            return self._reference_file
        
        raise FileError("Reference FASTA file not available for normalization")
    
    def cleanup(self):
        """Clean up any temporary files created during preparation."""
        logger.debug("Cleaning up FilePreparator resources")
        if hasattr(self, 'prepared_files'):
            self.prepared_files.clear()