# This code is based on PLastZ.py, originally created by Antoine Ho (AntoineHo)
# and obtained from GitHub (https://github.com/AntoineHo/PLastZ).
# It has been modified and integrated into ddPrimer for direct integration
# and improved compatibility.
#
# Original code licensed under Apache License, Version 2.0
# Original copyright information preserved as per license requirements.
# Modified for integration with ddPrimer by Jakob (2025)

"""
LastZ Runner module for ddPrimer pipeline.

This module provides functionality for running LastZ alignments in parallel,
handling sequence splitting, alignment, and result aggregation.
"""

import os
import errno
import shutil
import subprocess
import shlex
import glob
import logging
from multiprocessing import Pool
from Bio import SeqIO

from ..config import Config
from ..config.exceptions import AlignmentError


class LastZRunner:
    """
    Integrated parallel LastZ runner for sequence alignment.
    Handles splitting, parallel alignment, and result aggregation.
    """
    
    def __init__(self, config=None):
        """
        Initialize with configuration settings.
        
        Args:
            config: Configuration object (defaults to global Config)
        """
        self.config = config if config else Config
        self.logger = logging.getLogger("ddPrimer.helpers")
    
    def create_directory(self, path):
        """
        Create a directory if it doesn't exist.
        
        Args:
            path (str): Directory path to create
            
        Returns:
            str: Absolute path to the directory
            
        Raises:
            AlignmentError: If directory cannot be created
        """
        try:
            os.makedirs(path, exist_ok=True)
            return os.path.abspath(path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                self.logger.error(f"Cannot create directory: {path}")
                self.logger.debug(f"Error details: {str(e)}", exc_info=True)
                raise AlignmentError(f"Cannot create directory: {path}")
            else:
                return os.path.abspath(path)
    
    def run_command(self, cmd):
        """
        Execute a shell command with error handling.
        
        Args:
            cmd (str): Command to execute
            
        Returns:
            int: Return code from the command
        """
        try:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            stdout, _ = proc.communicate()
            
            if proc.returncode != 0:
                self.logger.warning(f"Command returned non-zero exit code: {proc.returncode}")
                self.logger.debug(f"Command: {cmd}")
                self.logger.debug(f"Output: {stdout.decode() if stdout else 'No output'}")
            
            return proc.returncode
        except Exception as e:
            self.logger.error(f"Error executing command: {cmd}")
            self.logger.debug(f"Error details: {str(e)}", exc_info=True)
            return -1
    
    def run_commands_parallel(self, commands, processes=None):
        """
        Run multiple commands in parallel.
        
        Args:
            commands (list): List of commands to execute
            processes (int): Number of parallel processes (default: from config)
            
        Raises:
            AlignmentError: If commands fail to execute
        """
        if not commands:
            self.logger.warning("No commands to execute")
            return
            
        if processes is None:
            processes = self.config.NUM_PROCESSES
        
        self.logger.debug(f"Running {len(commands)} commands in parallel using {processes} processes")
        
        try:
            with Pool(processes=processes) as pool:
                results = pool.map(self.run_command, commands)
                
            # Check for failures
            failed = sum(1 for r in results if r != 0)
            if failed > 0:
                self.logger.warning(f"{failed}/{len(commands)} commands failed")
            else:
                self.logger.debug(f"All {len(commands)} commands completed successfully")
                
        except Exception as e:
            self.logger.error(f"Error in parallel command execution")
            self.logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise AlignmentError(f"Parallel command execution failed: {str(e)}")
    
    def create_alignment_jobs(self, query_path, target_path, temp_dir, lastz_options=None):
        """
        Create extraction and alignment jobs.
        
        Args:
            query_path (str): Path to query FASTA file
            target_path (str): Path to target FASTA file
            temp_dir (str): Directory to store temporary files
            lastz_options (str): Additional options for LastZ
            
        Returns:
            tuple: Lists of extraction and alignment job commands
            
        Raises:
            AlignmentError: If job creation fails
        """
        extract_jobs = []
        align_jobs = []
        extracted_sequences = []
        pairs_done = []
        
        self.logger.debug(f"Creating alignment jobs for {query_path} vs {target_path}")
        
        try:
            # Extract unique sequence IDs from query and target
            query_contigs = []
            for qry in SeqIO.parse(query_path, "fasta"):
                query_contigs.append(qry.id)
                
            target_contigs = []
            for tgt in SeqIO.parse(target_path, "fasta"):
                target_contigs.append(tgt.id)
            
            self.logger.debug(f"Found {len(query_contigs)} query sequences and {len(target_contigs)} target sequences")
            
            # Create jobs for each sequence pair
            for query_id in query_contigs:
                if query_id not in extracted_sequences:
                    # Extract query sequence
                    out_fa_path = os.path.join(temp_dir, f"{query_id}.fa")
                    query_quoted = shlex.quote(query_path)
                    out_quoted = shlex.quote(out_fa_path)
                    extract_jobs.append(f"samtools faidx {query_quoted} {query_id} > {out_quoted}")
                    extracted_sequences.append(query_id)
                
                for target_id in target_contigs:
                    # Skip if we've already processed this pair (in either direction)
                    if [query_id, target_id] in pairs_done:
                        continue
                    
                    # Mark both directions as processed
                    pairs_done.append([query_id, target_id])
                    pairs_done.append([target_id, query_id])
                    
                    # Extract target sequence if not already done
                    if target_id not in extracted_sequences:
                        out_fa_path = os.path.join(temp_dir, f"{target_id}.fa")
                        target_quoted = shlex.quote(target_path)
                        out_quoted = shlex.quote(out_fa_path)
                        extract_jobs.append(f"samtools faidx {target_quoted} {target_id} > {out_quoted}")
                        extracted_sequences.append(target_id)
                    
                    # Create alignment job
                    tmp_output = os.path.join(temp_dir, f"{query_id}_V_{target_id}.TMP")
                    query_fa = os.path.join(temp_dir, f"{query_id}.fa")
                    target_fa = os.path.join(temp_dir, f"{target_id}.fa")
                    
                    # Quote paths for shell safety
                    query_fa_quoted = shlex.quote(query_fa)
                    target_fa_quoted = shlex.quote(target_fa)
                    tmp_output_quoted = shlex.quote(tmp_output)
                    
                    # Build LastZ command
                    if lastz_options:
                        cmd = f"lastz {query_fa_quoted} {target_fa_quoted} {lastz_options} > {tmp_output_quoted}"
                    else:
                        cmd = f"lastz {query_fa_quoted} {target_fa_quoted} > {tmp_output_quoted}"
                    
                    align_jobs.append(cmd)
            
            self.logger.debug(f"Created {len(extract_jobs)} extraction jobs and {len(align_jobs)} alignment jobs")
            return extract_jobs, align_jobs
            
        except Exception as e:
            self.logger.error(f"Error creating alignment jobs")
            self.logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise AlignmentError(f"Failed to create alignment jobs: {str(e)}")
    
    def run_parallel_alignment(self, ref_path, qry_path, output_dir, lastz_options=None, processes=None, keep_temp=False):
        """
        Run LastZ alignment between reference and query genomes in parallel.
        
        Args:
            ref_path (str): Path to reference FASTA file
            qry_path (str): Path to query FASTA file
            output_dir (str): Directory to store results
            lastz_options (str): Additional options for LastZ
            processes (int): Number of parallel processes
            keep_temp (bool): Whether to keep temporary files
            
        Returns:
            str: Path to the output MAF file
            
        Raises:
            AlignmentError: If alignment fails
        """
        # Set default processes if not provided
        if processes is None:
            processes = self.config.NUM_PROCESSES
        
        # Ensure paths are absolute
        ref_path = os.path.abspath(ref_path)
        qry_path = os.path.abspath(qry_path)
        
        # Extract base filenames without extensions for output naming
        ref_basename = os.path.splitext(os.path.basename(ref_path))[0]
        qry_basename = os.path.splitext(os.path.basename(qry_path))[0]
        
        # Create output directory and subdirectories
        input_dir = os.path.dirname(ref_path)
        alignments_dir = self.create_directory(os.path.join(input_dir, "Alignments"))
        temp_dir = self.create_directory(os.path.join(output_dir, "temp"))
        
        # Track .fai files created by samtools
        fai_files = []
        if ref_path != qry_path:
            fai_files = [f"{ref_path}.fai", f"{qry_path}.fai"]
        else:
            fai_files = [f"{ref_path}.fai"]
        
        try:
            self.logger.debug("\nPreparing LastZ alignment...")
            extract_jobs, align_jobs = self.create_alignment_jobs(ref_path, qry_path, temp_dir, lastz_options)
            
            self.logger.debug(f"Extracting sequences from reference and query genomes...")
            self.run_commands_parallel(extract_jobs, processes)
            
            self.logger.info("Running LastZ alignments...")
            self.run_commands_parallel(align_jobs, processes)
            
            # Combine all alignment results with new naming format
            output_file = os.path.join(alignments_dir, f"{ref_basename}_vs_{qry_basename}.maf")
            self.logger.debug(f"Combining alignment results to {output_file}...")
            
            # Use Python's file I/O capabilities for a more elegant solution
            with open(output_file, 'w') as outfile:
                # Get all TMP files in the temp directory
                temp_files = glob.glob(os.path.join(temp_dir, "*.TMP"))
                for temp_file in temp_files:
                    with open(temp_file, 'r') as infile:
                        outfile.write(infile.read())
            
            # Clean up temporary files
            if not keep_temp:
                self._clean_up_temp_files(temp_dir, fai_files)
            
            self.logger.debug(f"LastZ alignment completed successfully: {output_file}")
            return output_file
            
        except Exception as e:
            self.logger.error(f"Error in LastZ alignment")
            self.logger.debug(f"Error details: {str(e)}", exc_info=True)
            
            # Attempt to clean up even if there was an error
            if not keep_temp:
                try:
                    self._clean_up_temp_files(temp_dir, fai_files)
                except Exception as cleanup_error:
                    self.logger.warning(f"Error during cleanup after alignment failure: {str(cleanup_error)}")
            
            raise AlignmentError(f"LastZ alignment failed: {str(e)}")
    

    def _clean_up_temp_files(self, temp_dir, fai_files):
        """
        Clean up temporary files created during alignment.
        
        Args:
            temp_dir (str): Temporary directory to remove
            fai_files (list): List of .fai index files to remove
        """
        self.logger.debug("Cleaning up temporary files...")
        
        # Remove temp directory
        if os.path.exists(temp_dir):
            try:
                shutil.rmtree(temp_dir)
                self.logger.debug(f"Removed temporary directory: {temp_dir}")
            except Exception as e:
                self.logger.warning(f"Could not remove temporary directory {temp_dir}: {str(e)}")
        
        # Clean up .fai files
        for fai_file in fai_files:
            if os.path.exists(fai_file):
                try:
                    os.remove(fai_file)
                    self.logger.debug(f"Removed index file: {fai_file}")
                except Exception as e:
                    self.logger.warning(f"Could not remove index file {fai_file}: {str(e)}")