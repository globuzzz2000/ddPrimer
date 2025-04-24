#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BLAST Database Creator

A module for creating BLAST databases from FASTA files.
"""

import os
import subprocess
import logging
import traceback
import shlex
import tempfile
import shutil

class BlastDBCreator:
    """Creates and manages BLAST databases from FASTA files."""
    
    def create_database(self, fasta_file, db_name=None, output_dir=None, logger=None):
        """
        Compatibility method to match the call in pipeline.py.
        
        Args:
            fasta_file (str): Path to the FASTA file
            db_name (str, optional): Custom name for the database
            output_dir (str, optional): Directory to store the database files
            logger (logging.Logger, optional): Logger instance to use
            
        Returns:
            str: Path to the created BLAST database
        """
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
            raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
            
        # Create output directory if needed
        if output_dir is None:
            output_dir = os.path.join(os.path.dirname(fasta_file), "blast_db")
        
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
            
            # Verify the database was created properly
            if not BlastDBCreator.verify_db(temp_db_path, logger):
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
            logger.debug(traceback.format_exc())
            
            # Clean up temp directory if it still exists
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
                
            raise
    
    @staticmethod
    def verify_db(db_path, logger=None):
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
            quoted_db_path = f'"{db_path}"'
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