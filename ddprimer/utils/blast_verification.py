#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BLAST Database Verification Utilities

This module provides functions to verify BLAST databases before pipeline execution.
"""

import os
import subprocess
import tempfile
import logging
import re
import shlex
import time
from ..config import Config
from .blast_db_creator import BlastDBCreator

def verify_blast_database(logger=None):
    """
    Verify that a valid BLAST database exists and can be accessed.
    
    Args:
        logger: Logger instance for output messages
        
    Returns:
        bool: True if database is valid, False otherwise
    """
    if logger is None:
        logger = logging.getLogger("ddPrimer")
    
    db_path = Config.DB_PATH
    
    # Ensure path is absolute and properly expanded
    db_path = os.path.abspath(os.path.expanduser(db_path))
    
    logger.debug(f"Checking BLAST database at: {db_path}")
    
    # Basic file check - do the database files exist?
    required_extensions = ['.nhr', '.nin', '.nsq']
    missing_files = []
    
    for ext in required_extensions:
        if not os.path.exists(db_path + ext):
            missing_files.append(ext)
    
    if missing_files:
        logger.warning(f"BLAST database at {db_path} is missing files: {', '.join(missing_files)}")
        # If basic files are missing, return False right away
        logger.error("\n======================================")
        logger.error("ERROR: BLAST database files missing!")
        logger.error("\nPossible solutions:")
        logger.error("1. Create a BLAST database with the makeblastdb command:")
        logger.error("   makeblastdb -in your_genome.fasta -dbtype nucl -out your_db_name")
        logger.error("\n2. Create a BLAST database directly with ddprimer:")
        logger.error("   python -m ddprimer --dbfasta your_genome.fasta [--dbname custom_name] [--dboutdir output_dir]")
        logger.error("\n3. Set a different database path in your configuration file:")
        logger.error("   python -m ddprimer --config your_config.json")
        logger.error("======================================")
        return False
    else:
        logger.debug("BLAST database files found.")
    
    # Run a simple BLAST command to test functionality
    try:
        logger.debug("Running test BLAST command...")
        
        # Create a small temporary FASTA file
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp_query:
            tmp_query.write(">test_seq\nACGTACGTACGTACGTACGT\n")
            tmp_query.flush()
            tmp_filename = tmp_query.name
        
        try:
            # Quote the database path exactly like blast_processor.py does
            quoted_db_path = f'"{db_path}"'
            
            # Run test BLAST command with a very simple query
            cmd = [
                "blastn",
                "-task", "blastn-short",
                "-db", quoted_db_path,  # Use the quoted path here
                "-query", tmp_filename,
                "-outfmt", "6",
                "-max_target_seqs", "5"  # Increased to avoid warning
            ]
            
            # Log the command with proper quoting for debug purposes
            cmd_str = " ".join(shlex.quote(str(c)) for c in cmd)
            logger.debug(f"Running test command: {cmd_str}")
            
            result = subprocess.run(
                cmd,
                text=True,
                capture_output=True,
                timeout=15  # Increased timeout
            )
            
            # Clean up temp file
            os.remove(tmp_filename)
            
            # Analyze the result
            if result.returncode == 0:
                # Even with no hits, a successful command should not have an error exit code
                logger.debug("BLAST database test command executed successfully")
                
                # Check if there were hits or just no matches in the database
                if result.stdout.strip():
                    logger.debug("Test BLAST command found hits in the database")
                else:
                    logger.debug("Test BLAST command executed but found no hits (this is normal for a test query)")
                    
                return True
            else:
                # Check for common errors that might still allow the database to work
                error_msg = result.stderr.strip()
                
                # Memory map errors can sometimes be ignored if the test command runs
                if "memory map file error" in error_msg:
                    logger.warning(f"BLAST database memory map warning: {error_msg}")
                    
                    # Check for warnings about examining matches which aren't fatal
                    if re.search(r"Examining \d+ or more matches is recommended", error_msg):
                        logger.warning("Non-critical BLAST warning detected, database may still work")
                        return True
                    
                    # Additional check - try again after a small delay
                    # Sometimes memory mapping issues resolve themselves
                    logger.warning("Memory map error detected, retrying after short delay...")
                    time.sleep(2)
                    
                    # Try one more time with a different approach
                    return _retry_verify_with_blastdbcmd(db_path, logger)
                    
                # Log the full error
                logger.error(f"BLAST database test failed: {error_msg}")
                return False
                
        except subprocess.TimeoutExpired:
            logger.error("BLAST command timed out - database may be corrupted or inaccessible")
            # Clean up temp file
            os.remove(tmp_filename)
            return False
            
    except Exception as e:
        logger.error(f"Error testing BLAST database: {str(e)}")
        return False
        
    # This code should not be reached, but just in case
    return _retry_verify_with_blastdbcmd(db_path, logger)

def _retry_verify_with_blastdbcmd(db_path, logger):
    """
    Fallback verification using blastdbcmd instead of blastn.
    
    Args:
        db_path (str): Path to the BLAST database
        logger: Logger instance
        
    Returns:
        bool: True if verification succeeds, False otherwise
    """
    if logger is None:
        logger = logging.getLogger("ddPrimer")
        
    try:
        # Quote the database path exactly like blast_processor.py does
        quoted_db_path = f'"{db_path}"'
        
        # Try running blastdbcmd instead of blastn
        logger.debug("Attempting fallback verification with blastdbcmd...")
        
        cmd = [
            "blastdbcmd",
            "-db", quoted_db_path,  # Use the quoted path here
            "-info"
        ]
        
        # Log the command with proper quoting
        cmd_str = " ".join(shlex.quote(str(c)) for c in cmd)
        logger.debug(f"Running fallback command: {cmd_str}")
        
        result = subprocess.run(
            cmd,
            text=True,
            capture_output=True,
            timeout=8
        )
        
        if result.returncode == 0:
            logger.debug("BLAST database verified with blastdbcmd")
            return True
        else:
            logger.error(f"blastdbcmd verification failed: {result.stderr}")
            
            # If we've exhausted all options
            logger.error("\n======================================")
            logger.error("ERROR: BLAST database verification failed!")
            
            if "memory map file error" in result.stderr:
                logger.error("Memory map error detected. This can happen due to:")
                logger.error("1. Paths with spaces in them")
                logger.error("2. Permission issues")
                logger.error("3. System resource limitations")
                logger.error("\nRecommended solutions:")
                logger.error("- Create the database in a path without spaces:")
                logger.error("  mkdir -p ~/blast_dbs")
                logger.error("  makeblastdb -in your_genome.fasta -dbtype nucl -out ~/blast_dbs/your_db_name")
                logger.error("- Update your config to use this new database:")
                logger.error("  DB_PATH: ~/blast_dbs/your_db_name")
            else:
                logger.error("The database exists but appears to be corrupted or inaccessible.")
                logger.error("\nTry rebuilding the database with:")
                logger.error("  makeblastdb -in your_genome.fasta -dbtype nucl -out your_db_name")
                
            logger.error("======================================")
            
            return False
            
    except Exception as e:
        logger.error(f"Error during fallback verification: {str(e)}")
        return False