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
from multiprocessing import Pool
from Bio import SeqIO
from ..config import Config

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
    
    def create_directory(self, path):
        """
        Create a directory if it doesn't exist.
        
        Args:
            path (str): Directory path to create
            
        Returns:
            str: Absolute path to the directory
            
        Raises:
            Exception: If directory cannot be created
        """
        try:
            os.makedirs(path)
            return os.path.abspath(path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise Exception(f"ERROR: Cannot create directory: {path}")
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
            proc.communicate()
            return proc.returncode
        except Exception as e:
            print(f"Error executing command: {cmd}\n{e}")
            return -1
    
    def run_commands_parallel(self, commands, processes=None):
        """
        Run multiple commands in parallel.
        
        Args:
            commands (list): List of commands to execute
            processes (int): Number of parallel processes (default: from config)
        """
        if processes is None:
            processes = self.config.NUM_PROCESSES
        
        with Pool(processes=processes) as pool:
            pool.map(self.run_command, commands)
    
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
        """
        extract_jobs = []
        align_jobs = []
        extracted_sequences = []
        pairs_done = []
        
        self.config.debug(f"Creating alignment jobs for {query_path} vs {target_path}")
        
        # Extract unique sequence IDs from query and target
        query_contigs = []
        for qry in SeqIO.parse(query_path, "fasta"):
            query_contigs.append(qry.id)
            
        target_contigs = []
        for tgt in SeqIO.parse(target_path, "fasta"):
            target_contigs.append(tgt.id)
        
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
        
        self.config.debug(f"Created {len(extract_jobs)} extraction jobs and {len(align_jobs)} alignment jobs")
        return extract_jobs, align_jobs
    
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
        
        # Create output directory and Alignments directly inside it
        output_dir = self.create_directory(output_dir)
        alignments_dir = self.create_directory(os.path.join(output_dir, "Alignments"))
        temp_dir = self.create_directory(os.path.join(output_dir, "temp"))
        
        # Track .fai files created by samtools
        fai_files = []
        if ref_path != qry_path:
            fai_files = [f"{ref_path}.fai", f"{qry_path}.fai"]
        else:
            fai_files = [f"{ref_path}.fai"]
        
        self.config.debug("\nPreparing LastZ alignment...")
        extract_jobs, align_jobs = self.create_alignment_jobs(ref_path, qry_path, temp_dir, lastz_options)
        
        self.config.debug(f"Extracting sequences from reference and query genomes...")
        self.run_commands_parallel(extract_jobs, processes)
        
        print("Running LastZ alignments...")
        self.run_commands_parallel(align_jobs, processes)
        
        # Combine all alignment results with new naming format
        output_file = os.path.join(alignments_dir, f"{ref_basename}_vs_{qry_basename}.maf")
        self.config.debug(f"Combining alignment results to {output_file}...")
        
        # Use Python's file I/O capabilities for a more elegant solution
        with open(output_file, 'w') as outfile:
            # Get all TMP files in the temp directory
            temp_files = glob.glob(os.path.join(temp_dir, "*.TMP"))
            for temp_file in temp_files:
                with open(temp_file, 'r') as infile:
                    outfile.write(infile.read())
        
        # Clean up temporary files
        if not keep_temp:
            self.config.debug("Cleaning up temporary files...")
            shutil.rmtree(temp_dir)
            
            # Clean up .fai files
            for fai_file in fai_files:
                if os.path.exists(fai_file):
                    try:
                        os.remove(fai_file)
                    except OSError as e:
                        print(f"Warning: Could not remove index file {fai_file}: {e}")
        
        return output_file


# Example usage if run as standalone script
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Run LastZ alignments in parallel")
    parser.add_argument("--ref", required=True, help="Reference FASTA file")
    parser.add_argument("--qry", required=True, help="Query FASTA file")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--options", help="LastZ options", default="")
    parser.add_argument("--processes", type=int, help="Number of parallel processes")
    parser.add_argument("--keep-temp", action="store_true", help="Keep temporary files")
    parser.add_argument("--debug", action="store_true", help="Enable debug mode")
    
    args = parser.parse_args()
    
    # Set debug mode if requested
    if args.debug:
        Config.DEBUG_MODE = True
    
    # Create runner and execute alignment
    runner = LastZRunner()
    output_file = runner.run_parallel_alignment(
        args.ref, 
        args.qry, 
        args.out, 
        args.options, 
        args.processes, 
        args.keep_temp
    )
    
    print(f"Alignment completed. Output file: {output_file}")