#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BLAST Database Verification Utilities

This module provides functions to verify BLAST databases before pipeline execution
and offers options to create databases from model organisms when verification fails.
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
from .model_organism_manager import ModelOrganismManager

class BlastVerification:
    """Handles BLAST database verification before pipeline execution."""
    
    @staticmethod
    def verify_blast_database(logger=None):
        """
        Verify that a valid BLAST database exists and can be accessed.
        If verification fails, offer options to create a new database.
        
        Args:
            logger (logging.Logger, optional): Logger instance for output messages
            
        Returns:
            bool: True if database is valid, False otherwise
        """
        if logger is None:
            logger = logging.getLogger("ddPrimer")
        
        # Initialize the Config singleton to ensure DB_PATH is properly set
        config = Config.get_instance()
        db_path = config.DB_PATH
        
        # If no database path is set, immediately offer to create one
        if db_path is None:
            logger.info("No BLAST database configured. You need to create or select one.")
            return BlastVerification._handle_failed_verification(logger)
        
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
            # Database files are missing, offer to create a new database
            return BlastVerification._handle_failed_verification(logger)
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
                # Run test BLAST command with a very simple query
                cmd = [
                    "blastn",
                    "-task", "blastn-short",
                    "-db", f'"{db_path}"',  # Use the quoted path here
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
                        return BlastVerification._retry_verify_with_blastdbcmd(db_path, logger)
                        
                    # Log the full error
                    logger.error(f"BLAST database test failed: {error_msg}")
                    # Database test failed, offer to create a new database
                    return BlastVerification._handle_failed_verification(logger)
                    
            except subprocess.TimeoutExpired:
                logger.error("BLAST command timed out - database may be corrupted or inaccessible")
                # Clean up temp file
                os.remove(tmp_filename)
                # Database test failed, offer to create a new database
                return BlastVerification._handle_failed_verification(logger)
                
        except Exception as e:
            logger.error(f"Error testing BLAST database: {str(e)}")
            logger.debug(str(e), exc_info=True)
            # Database test failed, offer to create a new database
            return BlastVerification._handle_failed_verification(logger)
            
        # This code should not be reached, but just in case
        return BlastVerification._retry_verify_with_blastdbcmd(db_path, logger)

    @staticmethod
    def _retry_verify_with_blastdbcmd(db_path, logger):
        """
        Fallback verification using blastdbcmd instead of blastn.
        
        Args:
            db_path (str): Path to the BLAST database
            logger (logging.Logger): Logger instance
            
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
                # Database verification failed with all methods, offer to create a new database
                return BlastVerification._handle_failed_verification(logger)
                
        except Exception as e:
            logger.error(f"Error during fallback verification: {str(e)}")
            logger.debug(str(e), exc_info=True)
            # Database verification failed with all methods, offer to create a new database
            return BlastVerification._handle_failed_verification(logger)
    
    @staticmethod
    def _handle_failed_verification(logger):
        """
        Handle failed database verification by offering options:
        1. Create a new database from a model organism
        2. Create a new database from a custom file
        3. Try again with the current database
        4. Exit
        
        Args:
            logger (logging.Logger): Logger instance
            
        Returns:
            bool: True if a new database was created and verified, False otherwise
        """
        logger.error("\n======================================")
        logger.error("ERROR: BLAST database verification failed!")
        logger.error("======================================")
        logger.info("You have the following options:")
        logger.info("1. Create a new database from a model organism genome")
        logger.info("2. Create a new database from a custom FASTA file")
        logger.info("3. Try again with the current database")
        logger.info("4. Exit")
        
        try:
            choice = input("Enter your choice [1-4]: ")
            
            if choice == "1":
                # Create a new database from a model organism
                logger.info("Creating a new database from a model organism genome...")
                
                # Use ModelOrganismManager to select and fetch a model organism
                from .model_organism_manager import ModelOrganismManager
                organism_key, organism_name, fasta_file = ModelOrganismManager.select_model_organism(logger)
                
                if fasta_file is None:
                    logger.info("Database creation canceled. Exiting...")
                    return False
                
                # Use database name based on organism name
                if organism_key is not None:
                    organism_name = ModelOrganismManager.MODEL_ORGANISMS[organism_key]["name"]
                    # Extract just the scientific name (remove parentheses part)
                    scientific_name = organism_name.split(' (')[0] if ' (' in organism_name else organism_name
                    # Replace spaces with underscores
                    db_name = scientific_name.replace(' ', '_')
                    logger.debug(f"Using database name: {db_name}")
                else:
                    db_name = None
                
                # Create the database
                from .blast_db_creator import BlastDBCreator
                blast_db_creator = BlastDBCreator()
                db_path = blast_db_creator.create_database(fasta_file, db_name)
                
                # If a database was successfully created
                if db_path:
                    # Get the config instance to update the database path
                    config = Config.get_instance()
                    
                    # Check if there's already a database path set that's not the default
                    if config.DB_PATH and config.DB_PATH != "/Library/Application Support/Blast_DBs/Tair DB/TAIR10":
                        # If there's already a non-default database, ask if we should use the new one
                        logger.info(f"Current BLAST database path: {config.DB_PATH}")
                        use_new_db = input("Use the newly created database instead? [Y/n]: ").strip().lower()
                        
                        if use_new_db == "" or use_new_db.startswith("y"):
                            config.DB_PATH = db_path
                            config.USE_CUSTOM_DB = True
                            # Save the database path for future runs
                            Config.save_database_config(db_path)
                            logger.info(f"Now using new BLAST database: {db_path}")
                        else:
                            logger.info(f"Keeping current BLAST database: {config.DB_PATH}")
                    else:
                        # If no database or default database, automatically use the new one
                        config.DB_PATH = db_path
                        config.USE_CUSTOM_DB = True
                        # Save the database path for future runs
                        Config.save_database_config(db_path)
                        logger.debug(f"BLAST database created and set as active: {db_path}")
                    
                    # Clean up genome file if it was from a model organism
                    if fasta_file and organism_key is not None:
                        ModelOrganismManager.cleanup_genome_file(fasta_file, logger)
                    
                    return True
                else:
                    logger.error("Failed to create database. Please try again.")
                    return False
                    
            elif choice == "2":
                # Create a new database from a custom FASTA file
                logger.info("Creating a new database from a custom FASTA file...")
                
                # Create the database using BlastDBCreator
                from .blast_db_creator import BlastDBCreator
                blast_db_creator = BlastDBCreator()
                db_path = blast_db_creator.create_database(True)  # True means prompt for file
                
                if db_path is None:
                    logger.info("Database creation canceled. Exiting...")
                    return False
                
                # If a database was successfully created
                if db_path:
                    # Get the config instance to update the database path
                    config = Config.get_instance()
                    
                    # Check if there's already a database path set that's not the default
                    if config.DB_PATH and config.DB_PATH != "/Library/Application Support/Blast_DBs/Tair DB/TAIR10":
                        # If there's already a non-default database, ask if we should use the new one
                        logger.info(f"Current BLAST database path: {config.DB_PATH}")
                        use_new_db = input("Use the newly created database instead? [Y/n]: ").strip().lower()
                        
                        if use_new_db == "" or use_new_db.startswith("y"):
                            config.DB_PATH = db_path
                            config.USE_CUSTOM_DB = True
                            # Save the database path for future runs
                            Config.save_database_config(db_path)
                            logger.info(f"Now using new BLAST database: {db_path}")
                        else:
                            logger.info(f"Keeping current BLAST database: {config.DB_PATH}")
                    else:
                        # If no database or default database, automatically use the new one
                        config.DB_PATH = db_path
                        config.USE_CUSTOM_DB = True
                        # Save the database path for future runs
                        Config.save_database_config(db_path)
                        logger.info(f"BLAST database created and set as active: {db_path}")
                    
                    return True
                else:
                    logger.error("Failed to create database. Please try again.")
                    return False
                    
            elif choice == "3":
                # Try again with the current database
                logger.info("Retrying with the current database...")
                logger.info("Please ensure all BLAST processes are closed and try again.")
                return False
                
            elif choice == "4" or choice.lower() == "exit" or choice.lower() == "quit":
                # Exit
                logger.info("Exiting...")
                return False
                
            else:
                logger.error("Invalid choice. Please enter a number between 1 and 4.")
                # For invalid choices, recursively call this method
                # This line is important - when overriding in tests, we need to ensure 
                # we use the class name explicitly to maintain the patching chain
                return BlastVerification._handle_failed_verification(logger)
                
        except KeyboardInterrupt:
            logger.info("\nOperation canceled by user.")
            return False
        except Exception as e:
            logger.error(f"Error handling failed verification: {str(e)}")
            logger.debug(str(e), exc_info=True)
            return False