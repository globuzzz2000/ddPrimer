#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BLAST Database Verification Utilities for ddPrimer pipeline.

Contains functionality for:
1. BLAST database verification and validation
2. Database file integrity checking
3. Model organism database creation workflow
4. Fallback verification methods

This module provides functions to verify BLAST databases before pipeline
execution and offers options to create databases from model organisms
when verification fails.
"""

import os
import subprocess
import tempfile
import re
import shlex
import time
import logging
from ..config import Config, ExternalToolError, FileError

# Set up module logger
logger = logging.getLogger(__name__)


class BlastVerification:
    """
    Handles BLAST database verification before pipeline execution.
    
    This class provides comprehensive database verification methods
    including file integrity checks, functional testing, and fallback
    verification strategies when primary methods fail.
    
    Example:
        >>> if BlastVerification.verify_blast_database():
        ...     print("Database verified successfully")
        ... else:
        ...     print("Database verification failed")
    """
    
    @staticmethod
    def verify_blast_database():
        """
        Verify that a valid BLAST database exists and can be accessed.
        
        Performs comprehensive verification including file existence checks
        and functional testing. If verification fails, offers options to
        create a new database from model organisms.
        
        Returns:
            True if database is valid, False otherwise
            
        Raises:
            ExternalToolError: If BLAST tools are not available or functional
            FileError: If database files cannot be accessed
        """
        logger.debug("=== BLAST DATABASE VERIFICATION DEBUG ===")
        
        config = Config.get_instance()
        db_path = Config.DB_PATH
        logger.debug(f"Checking database path: {db_path}")
        
        # If no database path is set, immediately offer to create one
        if db_path is None:
            logger.info("No BLAST database configured. You need to create or select one.")
            logger.debug("=== END BLAST DATABASE VERIFICATION DEBUG ===")
            return BlastVerification._handle_failed_verification()
        
        # Ensure path is absolute and properly expanded
        db_path = os.path.abspath(os.path.expanduser(db_path))
        logger.debug(f"Expanded database path: {db_path}")
        
        # Basic file check - do the database files exist?
        required_extensions = ['.nhr', '.nin', '.nsq']
        missing_files = []
        
        for ext in required_extensions:
            file_path = db_path + ext
            if not os.path.exists(file_path):
                missing_files.append(ext)
            else:
                logger.debug(f"Found database file: {file_path}")
        
        if missing_files:
            logger.warning(f"BLAST database at {db_path} is missing files: {', '.join(missing_files)}")
            logger.debug("=== END BLAST DATABASE VERIFICATION DEBUG ===")
            return BlastVerification._handle_failed_verification()
        else:
            logger.debug("All required BLAST database files found")
        
        # Run a simple BLAST command to test functionality
        try:
            logger.debug("Running test BLAST command...")
            
            # Create a small temporary FASTA file
            try:
                with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp_query:
                    tmp_query.write(">test_seq\nACGTACGTACGTACGTACGT\n")
                    tmp_query.flush()
                    tmp_filename = tmp_query.name
            except OSError as e:
                error_msg = f"Failed to create temporary test file: {str(e)}"
                logger.error(error_msg)
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                raise FileError(error_msg) from e
            
            try:
                # Run test BLAST command with a very simple query
                cmd = [
                    "blastn",
                    "-task", "blastn-short",
                    "-db", f'"{db_path}"',
                    "-query", tmp_filename,
                    "-outfmt", "6",
                    "-max_target_seqs", "5"
                ]
                
                # Log the command with proper quoting for debug purposes
                cmd_str = " ".join(shlex.quote(str(c)) for c in cmd)
                logger.debug(f"Running test command: {cmd_str}")
                
                result = subprocess.run(
                    cmd,
                    text=True,
                    capture_output=True,
                    timeout=15
                )
                
                # Clean up temp file
                os.remove(tmp_filename)
                
                # Analyze the result
                if result.returncode == 0:
                    logger.debug("BLAST database test command executed successfully")
                    
                    # Check if there were hits or just no matches in the database
                    if result.stdout.strip():
                        logger.debug("Test BLAST command found hits in the database")
                    else:
                        logger.debug("Test BLAST command executed but found no hits (this is normal for a test query)")
                        
                    logger.debug("=== END BLAST DATABASE VERIFICATION DEBUG ===")
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
                            logger.debug("=== END BLAST DATABASE VERIFICATION DEBUG ===")
                            return True
                        
                        # Additional check - try again after a small delay
                        logger.warning("Memory map error detected, retrying after short delay...")
                        time.sleep(2)
                        
                        # Try one more time with a different approach
                        result = BlastVerification._retry_verify_with_blastdbcmd(db_path)
                        logger.debug("=== END BLAST DATABASE VERIFICATION DEBUG ===")
                        return result
                        
                    # Log the full error
                    error_msg = f"BLAST database test failed: {error_msg}"
                    logger.error(error_msg)
                    logger.debug("=== END BLAST DATABASE VERIFICATION DEBUG ===")
                    raise ExternalToolError(error_msg, tool_name="blastn")
                    
            except subprocess.TimeoutExpired:
                error_msg = "BLAST command timed out - database may be corrupted or inaccessible"
                logger.error(error_msg)
                # Clean up temp file
                try:
                    os.remove(tmp_filename)
                except OSError:
                    pass
                logger.debug("=== END BLAST DATABASE VERIFICATION DEBUG ===")
                raise ExternalToolError(error_msg, tool_name="blastn")
                
        except (ExternalToolError, FileError):
            # Re-raise custom exceptions without modification
            raise
        except Exception as e:
            error_msg = f"Error testing BLAST database: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END BLAST DATABASE VERIFICATION DEBUG ===")
            raise ExternalToolError(error_msg, tool_name="blastn") from e
            
        # This code should not be reached, but just in case
        result = BlastVerification._retry_verify_with_blastdbcmd(db_path)
        logger.debug("=== END BLAST DATABASE VERIFICATION DEBUG ===")
        return result

    @staticmethod
    def _retry_verify_with_blastdbcmd(db_path):
        """
        Fallback verification using blastdbcmd instead of blastn.
        
        Provides an alternative verification method when standard BLAST
        testing fails, using database information commands instead of
        actual sequence searches.
        
        Args:
            db_path: Path to the BLAST database
            
        Returns:
            True if verification succeeds, False otherwise
            
        Raises:
            ExternalToolError: If blastdbcmd is not available or fails
        """
        try:
            # Quote the database path exactly like blast_processor.py does
            quoted_db_path = f'"{db_path}"'
            
            # Try running blastdbcmd instead of blastn
            logger.debug("Attempting fallback verification with blastdbcmd...")
            
            cmd = [
                "blastdbcmd",
                "-db", quoted_db_path,
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
                error_msg = f"blastdbcmd verification failed"
                logger.error(error_msg)
                logger.debug(f"blastdbcmd stderr: {result.stderr}")
                raise ExternalToolError(error_msg, tool_name="blastdbcmd")
                
        except subprocess.TimeoutExpired:
            error_msg = "blastdbcmd verification timed out"
            logger.error(error_msg)
            raise ExternalToolError(error_msg, tool_name="blastdbcmd")
        except ExternalToolError:
            # Re-raise ExternalToolError without modification
            raise
        except Exception as e:
            error_msg = f"Error during fallback verification: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ExternalToolError(error_msg, tool_name="blastdbcmd") from e
    
    @staticmethod
    def _handle_failed_verification():
        """
        Handle failed database verification by offering creation options.
        
        Provides interactive menu for users to:
        1. Create database from model organism
        2. Create database from custom file  
        3. Retry with current database
        4. Exit pipeline
        
        Returns:
            True if new database was created and verified, False otherwise
            
        Raises:
            ExternalToolError: If database creation tools fail
            FileError: If file operations fail
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
                try:
                    from . import ModelOrganismManager
                    organism_key, organism_name, fasta_file = ModelOrganismManager.select_model_organism()
                except Exception as e:
                    error_msg = f"Failed to access model organism selection: {str(e)}"
                    logger.error(error_msg)
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    raise ExternalToolError(error_msg, tool_name="ModelOrganismManager") from e
                
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
                try:
                    from . import BlastDBCreator
                    blast_db_creator = BlastDBCreator()
                    db_path = blast_db_creator.create_database(fasta_file, db_name)
                except Exception as e:
                    error_msg = f"Failed to create database: {str(e)}"
                    logger.error(error_msg)
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    raise ExternalToolError(error_msg, tool_name="BlastDBCreator") from e
                
                # If a database was successfully created
                if db_path:
                    # Check if there's already a database path set that's not the default
                    if Config.DB_PATH and Config.DB_PATH != "/Library/Application Support/Blast_DBs/Tair DB/TAIR10":
                        logger.info(f"Current BLAST database path: {Config.DB_PATH}")
                        use_new_db = input("Use the newly created database instead? [Y/n]: ").strip().lower()
                        
                        if use_new_db == "" or use_new_db.startswith("y"):
                            Config.DB_PATH = db_path
                            Config.USE_CUSTOM_DB = True
                            Config.save_database_config(db_path)
                            logger.info(f"Now using new BLAST database: {db_path}")
                        else:
                            logger.info(f"Keeping current BLAST database: {Config.DB_PATH}")
                    else:
                        Config.DB_PATH = db_path
                        Config.USE_CUSTOM_DB = True
                        Config.save_database_config(db_path)
                        logger.debug(f"BLAST database created and set as active: {db_path}")
                    
                    # Clean up genome file if it was from a model organism
                    if fasta_file and organism_key is not None:
                        try:
                            ModelOrganismManager.cleanup_genome_file(fasta_file)
                        except Exception as e:
                            logger.warning(f"Failed to clean up genome file: {str(e)}")
                            logger.debug(f"Cleanup error details: {str(e)}", exc_info=True)
                    
                    return True
                else:
                    logger.error("Failed to create database. Please try again.")
                    return False
                    
            elif choice == "2":
                # Create a new database from a custom FASTA file
                logger.info("Creating a new database from a custom FASTA file...")
                
                # Create the database using BlastDBCreator
                try:
                    from . import BlastDBCreator
                    blast_db_creator = BlastDBCreator()
                    db_path = blast_db_creator.create_database(True)  # True means prompt for file
                except Exception as e:
                    error_msg = f"Failed to create database: {str(e)}"
                    logger.error(error_msg)
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    raise ExternalToolError(error_msg, tool_name="BlastDBCreator") from e
                
                if db_path is None:
                    logger.info("Database creation canceled. Exiting...")
                    return False
                
                # If a database was successfully created
                if db_path:
                    # Check if there's already a database path set that's not the default
                    if Config.DB_PATH and Config.DB_PATH != "/Library/Application Support/Blast_DBs/Tair DB/TAIR10":
                        logger.info(f"Current BLAST database path: {Config.DB_PATH}")
                        use_new_db = input("Use the newly created database instead? [Y/n]: ").strip().lower()
                        
                        if use_new_db == "" or use_new_db.startswith("y"):
                            Config.DB_PATH = db_path
                            Config.USE_CUSTOM_DB = True
                            Config.save_database_config(db_path)
                            logger.info(f"Now using new BLAST database: {db_path}")
                        else:
                            logger.info(f"Keeping current BLAST database: {Config.DB_PATH}")
                    else:
                        Config.DB_PATH = db_path
                        Config.USE_CUSTOM_DB = True
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
                return BlastVerification._handle_failed_verification()
                
        except KeyboardInterrupt:
            logger.info("\nOperation canceled by user.")
            return False
        except (ExternalToolError, FileError):
            # Re-raise custom exceptions without modification
            raise
        except Exception as e:
            error_msg = f"Error handling failed verification: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ExternalToolError(error_msg, tool_name="database_creation") from e