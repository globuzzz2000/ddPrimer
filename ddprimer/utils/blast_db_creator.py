#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BLAST Database Creator for ddPrimer pipeline.

Contains functionality for:
1. BLAST database creation from FASTA files
2. Model organism genome integration
3. Database verification and validation
4. Temporary file management

This module integrates with the broader ddPrimer pipeline to provide
robust BLAST database management capabilities.
"""

import os
import subprocess
import shlex
import tempfile
import shutil
import logging
from pathlib import Path
from .model_organism_manager import ModelOrganismManager
from ..config import FileError, ExternalToolError

# Set up module logger
logger = logging.getLogger(__name__)


class BlastDBCreator:
    """
    Creates and manages BLAST databases from FASTA files.
    
    This class provides methods for creating BLAST databases from various
    sources including local FASTA files and model organism genomes. It handles
    temporary file management and database verification.
    
    Example:
        >>> creator = BlastDBCreator()
        >>> db_path = creator.create_database("genome.fasta", "my_db")
        >>> print(f"Database created at: {db_path}")
    """
    
    @staticmethod
    def create_db(fasta_file, output_dir=None, db_name=None):
        """
        Create a BLAST database from a FASTA file.
        
        Executes makeblastdb command to create a nucleotide BLAST database
        from the provided FASTA file with proper error handling and verification.
        
        Args:
            fasta_file: Path to the FASTA file
            output_dir: Directory to store database files
            db_name: Custom name for the database
            
        Returns:
            Path to the created BLAST database
            
        Raises:
            FileError: If FASTA file doesn't exist or output directory issues
            ExternalToolError: If makeblastdb execution fails
        """
        logger.debug(f"Creating BLAST database from {fasta_file}...")
        
        # Handle path with spaces
        fasta_file = os.path.abspath(os.path.expanduser(fasta_file))
        
        # Verify the FASTA file exists
        if not os.path.exists(fasta_file):
            error_msg = f"FASTA file not found: {fasta_file}"
            logger.error(error_msg)
            raise FileError(error_msg)
            
        # Create output directory if needed
        if output_dir is None:
            # Use system directory if running as root, otherwise use user home
            if os.geteuid() == 0:
                output_dir = "/usr/local/share/ddprimer/blast_db"
            else:
                output_dir = os.path.join(Path.home(), ".ddprimer", "blast_db")
        
        # Make sure output_dir is absolute path
        output_dir = os.path.abspath(os.path.expanduser(output_dir))
        try:
            os.makedirs(output_dir, exist_ok=True)
        except OSError as e:
            error_msg = f"Failed to create output directory {output_dir}: {str(e)}"
            logger.error(error_msg)
            raise FileError(error_msg) from e
        
        # Create a temporary working directory with a safe path
        try:
            temp_dir = tempfile.mkdtemp(prefix="blastdb_")
            temp_fasta = os.path.join(temp_dir, os.path.basename(fasta_file))
            shutil.copy2(fasta_file, temp_fasta)
        except (OSError, IOError) as e:
            error_msg = f"Failed to create temporary working directory: {str(e)}"
            logger.error(error_msg)
            raise FileError(error_msg) from e

        # Generate safe DB name and path in the temp directory
        safe_db_name = db_name or os.path.splitext(os.path.basename(fasta_file))[0]
        safe_db_name = safe_db_name.replace(" ", "_").replace("-", "_")
        temp_db_path = os.path.join(temp_dir, safe_db_name)
        
        db_path = os.path.join(output_dir, safe_db_name)
        logger.debug(f"Database will be created at: {db_path}")
        
        try:
            # Run makeblastdb command with proper escaping of paths
            cmd = [
                "makeblastdb",
                "-in", temp_fasta,
                "-dbtype", "nucl",
                "-out", temp_db_path
            ]
            
            # Log the command but with quotes for human readability
            cmd_str = " ".join(shlex.quote(str(c)) for c in cmd)
            logger.debug(f"Running command: {cmd_str}")
            
            result = subprocess.run(
                cmd,
                text=True,
                capture_output=True
            )
            
            if result.returncode != 0:
                error_msg = f"BLAST database creation failed"
                logger.error(error_msg)
                logger.debug(f"makeblastdb stderr: {result.stderr}")
                logger.debug(f"makeblastdb stdout: {result.stdout}")
                raise ExternalToolError(error_msg, tool_name="makeblastdb")
                
            logger.info(f"\nBLAST database created successfully at: {db_path}")
            logger.info("")
            
            # Verify the database was created properly
            if not BlastDBCreator._verify_db(temp_db_path):
                error_msg = "BLAST database verification failed"
                logger.error(error_msg)
                raise ExternalToolError(error_msg, tool_name="makeblastdb")
                
            # Copy the resulting files to the target output location
            db_files = [f for f in os.listdir(temp_dir) if f.startswith(safe_db_name) and f != safe_db_name]
            for db_file in db_files:
                src = os.path.join(temp_dir, db_file)
                dst = os.path.join(output_dir, db_file)
                if os.path.exists(src):
                    shutil.copy2(src, dst)
                else:
                    logger.warning(f"Expected database file not found: {src}")

            # Clean up temp directory
            shutil.rmtree(temp_dir)
            
            return db_path
            
        except ExternalToolError:
            # Re-raise ExternalToolError without modification
            raise
        except Exception as e:
            error_msg = f"Failed to create BLAST database: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            
            # Clean up temp directory if it still exists
            if os.path.exists(temp_dir):
                try:
                    shutil.rmtree(temp_dir)
                except OSError as cleanup_error:
                    logger.warning(f"Error cleaning up temporary directory: {str(cleanup_error)}")
                
            raise ExternalToolError(error_msg, tool_name="makeblastdb") from e
    
    @staticmethod
    def _verify_db(db_path):
        """
        Verify that a BLAST database exists and is valid.
        
        Checks for required database files and validates the database
        using blastdbcmd to ensure it can be accessed properly.
        
        Args:
            db_path: Path to the BLAST database
            
        Returns:
            True if the database is valid, False otherwise
        """
        # Check if database files exist
        required_extensions = ['.nhr', '.nin', '.nsq']
        missing_files = []
        
        for ext in required_extensions:
            if not os.path.exists(db_path + ext):
                missing_files.append(ext)
        
        if missing_files:
            logger.warning(f"BLAST database at {db_path} is missing files: {', '.join(missing_files)}")
            return False
            
        # Run blastdbcmd to verify the database
        try:
            cmd = [
                "blastdbcmd",
                "-db", db_path,
                "-info"
            ]
            
            # Log the command with quotes for human readability
            cmd_str = " ".join(shlex.quote(str(c)) for c in cmd)
            logger.debug(f"Running verification command: {cmd_str}")
            
            result = subprocess.run(
                cmd,
                text=True,
                capture_output=True
            )
            
            if result.returncode != 0:
                logger.warning(f"BLAST database verification failed: {result.stderr}")
                return False
                
            logger.debug(f"BLAST database verified: {db_path}")
            return True
            
        except Exception as e:
            logger.warning(f"Error verifying BLAST database: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return False