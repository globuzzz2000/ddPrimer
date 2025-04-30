#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BLAST Database Creator

A module for creating BLAST databases from FASTA files,
with support for automatic fetching of model organism genomes.
"""

import os
import subprocess
import logging
import traceback
import shlex
import tempfile
import shutil
from pathlib import Path
from .model_organism_manager import ModelOrganismManager


class BlastDBCreator:
    """Creates and manages BLAST databases from FASTA files."""
    
    def create_database(self, fasta_file=None, db_name=None, output_dir=None, logger=None):
        """
        Create a BLAST database from a FASTA file or a model organism.
        
        Args:
            fasta_file (str, optional): Path to the FASTA file or None to prompt for selection
            db_name (str, optional): Custom name for the database
            output_dir (str, optional): Directory to store the database files
            logger (logging.Logger, optional): Logger instance to use
            
        Returns:
            str: Path to the created BLAST database
            
        Raises:
            FileNotFoundError: If the FASTA file doesn't exist
            Exception: If database creation fails
        """
        # Use default logger if none provided
        if logger is None:
            logger = logging.getLogger("ddPrimer")
            
        # If no fasta_file is provided, offer model organism selection
        if fasta_file is None or fasta_file is True:  # True indicates to prompt for selection
            logger.info("No FASTA file specified. You can select a model organism or provide a custom file.")
            
            organism_key, organism_name, fasta_file = ModelOrganismManager.select_model_organism(logger)
            
            # If selection was canceled or failed
            if fasta_file is None:
                logger.info("BLAST database creation canceled.")
                return None
                
            # If a model organism was selected, use its name as the default DB name
            if organism_key is not None and db_name is None:
                organism_name = ModelOrganismManager.MODEL_ORGANISMS[organism_key]["name"]
                scientific_name = organism_name.split(' (')[0] if ' (' in organism_name else organism_name
                db_name = scientific_name.replace(' ', '_')
                logger.info(f"Using database name: {db_name}")
        
        return self.create_db(fasta_file, output_dir, db_name, logger)
    
    @staticmethod
    def create_db(fasta_file, output_dir=None, db_name=None, logger=None):
        """
        Create a BLAST database from a FASTA file.
        
        Args:
            fasta_file (str): Path to the FASTA file
            output_dir (str, optional): Directory to store the database files
            db_name (str, optional): Custom name for the database (default: derived from filename)
            logger (logging.Logger, optional): Logger instance to use
            
        Returns:
            str: Path to the created BLAST database
        """
        # Use default logger if none provided
        if logger is None:
            logger = logging.getLogger("ddPrimer")
        
        logger.info(f"Creating BLAST database from {fasta_file}...")
        
        # Handle path with spaces
        fasta_file = os.path.abspath(os.path.expanduser(fasta_file))
        
        # Verify the FASTA file exists
        if not os.path.exists(fasta_file):
            logger.error(f"FASTA file not found: {fasta_file}")
            raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
            
        # Create output directory if needed
        if output_dir is None:
            # Use system directory if running as root, otherwise use user home
            if os.geteuid() == 0:  # Check if running as root
                output_dir = "/usr/local/share/ddprimer/blast_db"
            else:
                output_dir = os.path.join(Path.home(), ".ddprimer", "blast_db")
        
        # Make sure output_dir is absolute path
        output_dir = os.path.abspath(os.path.expanduser(output_dir))
        os.makedirs(output_dir, exist_ok=True)
        
        # Create a temporary working directory with a safe path
        temp_dir = tempfile.mkdtemp(prefix="blastdb_")
        temp_fasta = os.path.join(temp_dir, os.path.basename(fasta_file))
        shutil.copy2(fasta_file, temp_fasta)

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
                logger.error(f"Error creating BLAST database: {result.stderr}")
                logger.debug(f"Command output: {result.stdout}")
                raise Exception(f"BLAST database creation failed: {result.stderr}")
                
            logger.info(f"BLAST database created successfully at: {db_path}")
            logger.info("")
            
            # Verify the database was created properly
            if not BlastDBCreator._verify_db(temp_db_path, logger):
                logger.error("BLAST database verification failed")
                raise Exception("BLAST database verification failed")
                
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
            
        except Exception as e:
            logger.error(f"Failed to create BLAST database: {str(e)}")
            logger.debug(traceback.format_exc(), exc_info=True)
            
            # Clean up temp directory if it still exists
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
                
            raise
    
    @staticmethod
    def _verify_db(db_path, logger=None):
        """
        Verify that a BLAST database exists and is valid.
        
        Args:
            db_path (str): Path to the BLAST database
            logger (logging.Logger, optional): Logger instance to use
            
        Returns:
            bool: True if the database is valid, False otherwise
        """
        # Use default logger if none provided
        if logger is None:
            logger = logging.getLogger("ddPrimer")
            
        # Check if database files exist
        required_extensions = ['.nhr', '.nin', '.nsq']
        missing_files = []
        
        for ext in required_extensions:
            if not os.path.exists(db_path + ext):
                missing_files.append(ext)
        
        if missing_files:
            logger.warning(f"BLAST database at {db_path} is missing files: {', '.join(missing_files)}")
            return False
            
        # Run blastdbcmd to verify the database - using the same approach as blast_processor
        # by wrapping the path in quotes
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
            return False