#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unified BLAST Database Manager for ddPrimer pipeline.

Contains comprehensive functionality for:
1. BLAST database creation from FASTA files
2. Model organism genome fetching and processing
3. Database verification and validation
4. Interactive database selection menus
5. Persistent database configuration management

This module provides a complete solution for BLAST database management,
integrating all previously separate components into a cohesive system.
"""

import os
import sys
import subprocess
import shlex
import tempfile
import shutil
import logging
import urllib.request
import gzip
import time
import re
from pathlib import Path
from typing import Optional, Dict, Tuple, Union
from ..config import Config, FileError, ExternalToolError

# Set up module logger
logger = logging.getLogger(__name__)


class BlastDatabaseManager:
    """
    Unified manager for all BLAST database operations in ddPrimer pipeline.
    
    This class provides comprehensive BLAST database management including:
    - Creation from FASTA files and model organisms
    - Verification and validation
    - Interactive selection menus
    - Persistent configuration management
    
    Example:
        >>> manager = BlastDatabaseManager()
        >>> if manager.verify_database():
        ...     print("Database ready for use")
        ... else:
        ...     success = manager.setup_database_interactive()
    """
    
    # Model organisms with their genome download URLs
    MODEL_ORGANISMS = {
        "Thale cress": {
            "name": "Arabidopsis thaliana",
            "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz",
            "compressed": True
        },
        "Sand rock-cress": {
            "name": "Arabidopsis arenosa",
            "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/905/216/605/GCA_905216605.1_AARE701a/GCA_905216605.1_AARE701a_genomic.fna.gz",
            "compressed": True
        },
        "E. coli": {
            "name": "Escherichia coli K-12 MG1655",
            "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz",
            "compressed": True
        },
        "Baker's yeast": {
            "name": "Saccharomyces cerevisiae",
            "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz",
            "compressed": True
        },
        "Human": {
            "name": "Homo sapiens",
            "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz",
            "compressed": True
        },
        "Mouse": {
            "name": "Mus musculus",
            "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz",
            "compressed": True
        },
        "Rat": {
            "name": "Rattus norvegicus",
            "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/895/GCF_000001895.5_Rnor_6.0/GCF_000001895.5_Rnor_6.0_genomic.fna.gz",
            "compressed": True
        },
        "Fruit fly": {
            "name": "Drosophila melanogaster",
            "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz",
            "compressed": True
        },
        "Roundworm": {
            "name": "Caenorhabditis elegans",
            "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.fna.gz",
            "compressed": True
        }
    }
    
    def __init__(self):
        """Initialize the BLAST database manager."""
        self.logger = logger
        
    #############################################################################
    #                           PUBLIC API METHODS
    #############################################################################
    
    def verify_database(self) -> bool:
        """
        Verify that a valid BLAST database exists and can be accessed.
        
        Performs comprehensive verification including file existence checks
        and functional testing.
        
        Returns:
            True if database is valid, False otherwise
            
        Raises:
            ExternalToolError: If BLAST tools are not available or functional
            FileError: If database files cannot be accessed
        """
        logger.debug("=== BLAST DATABASE VERIFICATION ===")
        
        db_path = Config.DB_PATH
        logger.debug(f"Checking database path: {db_path}")
        
        # If no database path is set, verification fails
        if db_path is None:
            logger.info("No BLAST database configured.")
            return False
        
        # Ensure path is absolute and properly expanded
        db_path = os.path.abspath(os.path.expanduser(db_path))
        logger.debug(f"Expanded database path: {db_path}")
        
        # Check for required database files
        if not self._check_database_files(db_path):
            return False
        
        # Run functional test
        try:
            return self._test_database_functionality(db_path)
        except (ExternalToolError, FileError):
            raise
        except Exception as e:
            error_msg = f"Unexpected error during database verification: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ExternalToolError(error_msg, tool_name="blast_verification") from e
    
    def setup_database_interactive(self) -> bool:
        """
        Interactive database setup when verification fails.
        
        Provides menu-driven interface for:
        1. Creating database from model organism
        2. Creating database from custom file
        3. Selecting existing database
        4. Retrying current database
        
        Returns:
            True if database was successfully set up, False otherwise
        """
        logger.error("\n======================================")
        logger.error("ERROR: BLAST database verification failed!")
        logger.error("======================================")
        logger.info("You have the following options:")
        logger.info("1. Create database from model organism")
        logger.info("2. Create database from custom FASTA file")
        logger.info("3. Select existing database")
        logger.info("4. Retry with current database")
        logger.info("5. Exit")
        
        try:
            choice = input("\nEnter your choice [1-5]: ")
            
            if choice == "1":
                return self._create_from_model_organism()
            elif choice == "2":
                return self._create_from_custom_file()
            elif choice == "3":
                return self._select_existing_database()
            elif choice == "4":
                logger.info("Retrying with current database...")
                return False
            elif choice == "5":
                logger.info("Exiting...")
                return False
            else:
                logger.error("Invalid choice. Please enter a number between 1 and 5.")
                return self.setup_database_interactive()
                
        except KeyboardInterrupt:
            logger.info("\nOperation canceled by user.")
            return False
        except Exception as e:
            error_msg = f"Error during interactive setup: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return False
    
    def create_database(self, fasta_file: str, db_name: Optional[str] = None, 
                       output_dir: Optional[str] = None) -> str:
        """
        Create a BLAST database from a FASTA file.
        
        Args:
            fasta_file: Path to the FASTA file
            db_name: Custom name for the database
            output_dir: Directory to store database files
            
        Returns:
            Path to the created BLAST database
            
        Raises:
            FileError: If FASTA file doesn't exist or output directory issues
            ExternalToolError: If makeblastdb execution fails
        """
        logger.debug(f"Creating BLAST database from {fasta_file}")
        
        # Validate input file
        fasta_file = os.path.abspath(os.path.expanduser(fasta_file))
        if not os.path.exists(fasta_file):
            error_msg = f"FASTA file not found: {fasta_file}"
            logger.error(error_msg)
            raise FileError(error_msg)
        
        # Set up output directory
        if output_dir is None:
            output_dir = self._get_default_db_directory()
        
        output_dir = os.path.abspath(os.path.expanduser(output_dir))
        try:
            os.makedirs(output_dir, exist_ok=True)
        except OSError as e:
            error_msg = f"Failed to create output directory {output_dir}: {str(e)}"
            logger.error(error_msg)
            raise FileError(error_msg) from e
        
        # Generate database name and path
        safe_db_name = db_name or os.path.splitext(os.path.basename(fasta_file))[0]
        safe_db_name = safe_db_name.replace(" ", "_").replace("-", "_")
        db_path = os.path.join(output_dir, safe_db_name)
        
        # Create database using temporary directory for safety
        temp_dir = None
        try:
            temp_dir = tempfile.mkdtemp(prefix="blastdb_")
            temp_fasta = os.path.join(temp_dir, os.path.basename(fasta_file))
            shutil.copy2(fasta_file, temp_fasta)
            
            temp_db_path = os.path.join(temp_dir, safe_db_name)
            
            # Run makeblastdb
            cmd = [
                "makeblastdb",
                "-in", temp_fasta,
                "-dbtype", "nucl",
                "-out", temp_db_path
            ]
            
            cmd_str = " ".join(shlex.quote(str(c)) for c in cmd)
            logger.debug(f"Running command: {cmd_str}")
            
            result = subprocess.run(cmd, text=True, capture_output=True)
            
            if result.returncode != 0:
                error_msg = f"BLAST database creation failed"
                logger.error(error_msg)
                logger.debug(f"makeblastdb stderr: {result.stderr}")
                raise ExternalToolError(error_msg, tool_name="makeblastdb")
            
            # Verify the database was created
            if not self._check_database_files(temp_db_path):
                error_msg = "BLAST database verification failed after creation"
                logger.error(error_msg)
                raise ExternalToolError(error_msg, tool_name="makeblastdb")
            
            # Copy database files to final location
            db_files = [f for f in os.listdir(temp_dir) 
                       if f.startswith(safe_db_name) and f != safe_db_name]
            for db_file in db_files:
                src = os.path.join(temp_dir, db_file)
                dst = os.path.join(output_dir, db_file)
                if os.path.exists(src):
                    shutil.copy2(src, dst)
            
            logger.info(f"\nNew BLAST database created at: {db_path}")
            return db_path
            
        except (FileError, ExternalToolError):
            raise
        except Exception as e:
            error_msg = f"Failed to create BLAST database: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ExternalToolError(error_msg, tool_name="makeblastdb") from e
        finally:
            if temp_dir and os.path.exists(temp_dir):
                try:
                    shutil.rmtree(temp_dir)
                except OSError as cleanup_error:
                    logger.warning(f"Error cleaning up temporary directory: {str(cleanup_error)}")
    
    def select_model_organism(self) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        """
        Present menu to select a model organism or existing database.
        
        Returns:
            Tuple of (organism_key, organism_name, file_path) or special values:
            - For model organism: (organism_key, organism_name, fasta_file_path)
            - For custom file: (None, "Custom file", fasta_file_path)
            - For existing database: ('existing_db', database_name, database_path)
            - For cancellation: (None, None, None)
        """
        menu = self._get_organism_menu()
        logger.info(menu)
        
        try:
            choice = input(f"Enter your choice [0-{len(self.MODEL_ORGANISMS) + 2}]: ")
            choice = int(choice)
            
            if choice == 0:  # Custom file
                return self._handle_custom_file_selection()
            elif choice == len(self.MODEL_ORGANISMS) + 1:  # Existing database
                return self._handle_existing_database_selection()
            elif choice == len(self.MODEL_ORGANISMS) + 2:  # Cancel
                return None, None, None
            elif 1 <= choice <= len(self.MODEL_ORGANISMS):
                return self._handle_model_organism_selection(choice)
            else:
                logger.error("Invalid choice. Please enter a number within the range.")
                return None, None, None
                
        except (ValueError, KeyboardInterrupt):
            logger.info("Operation canceled.")
            return None, None, None
        except Exception as e:
            error_msg = f"Error in organism selection: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return None, None, None
    
    def find_existing_databases(self) -> Dict[str, str]:
        """
        Find all existing BLAST databases in standard locations.
        
        Returns:
            Dictionary mapping database paths to display names
        """
        logger.debug("Searching for existing BLAST databases")
        databases = {}
        
        # Standard locations to search
        search_locations = [
            "/usr/local/share/ddprimer/blast_db",
            os.path.join(Path.home(), ".ddprimer", "blast_db"),
            "/Library/Application Support/Blast_DBs"  # macOS location
        ]
        
        for location in search_locations:
            if os.path.exists(location):
                try:
                    for root, dirs, files in os.walk(location):
                        db_files = [f for f in files if f.endswith(".nhr")]
                        
                        for db_file in db_files:
                            db_name = os.path.splitext(db_file)[0]
                            db_path = os.path.join(root, db_name)
                            
                            # Create display name
                            rel_path = os.path.relpath(db_path, location)
                            display_name = rel_path.replace(os.sep, " / ").replace("_", " ")
                            
                            databases[db_path] = display_name
                            logger.debug(f"Found database: {display_name}")
                        
                except OSError as e:
                    logger.warning(f"Error searching {location}: {str(e)}")
        
        logger.debug(f"Found {len(databases)} existing databases")
        return databases
    
    def set_active_database(self, db_path: str) -> bool:
        """
        Set the active BLAST database and save configuration.
        
        Args:
            db_path: Path to the BLAST database
            
        Returns:
            True if successfully set, False otherwise
        """
        try:
            # Verify the database exists and is valid
            if not self._check_database_files(db_path):
                logger.error(f"Invalid database path: {db_path}")
                return False
            
            # Update configuration
            Config.DB_PATH = db_path
            Config.USE_CUSTOM_DB = True
            Config.save_database_config(db_path)
            
            logger.info(f"\nSelected BLAST database: {os.path.basename(db_path)}")
            return True
            
        except Exception as e:
            error_msg = f"Failed to set active database: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return False
    
    #############################################################################
    #                           PRIVATE HELPER METHODS
    #############################################################################
    
    def _check_database_files(self, db_path: str) -> bool:
        """Check if required BLAST database files exist."""
        required_extensions = ['.nhr', '.nin', '.nsq']
        missing_files = []
        
        for ext in required_extensions:
            if not os.path.exists(db_path + ext):
                missing_files.append(ext)
        
        if missing_files:
            logger.warning(f"BLAST database missing files: {', '.join(missing_files)}")
            return False
        
        logger.debug("All required BLAST database files found")
        return True
    
    def _test_database_functionality(self, db_path: str) -> bool:
        """Test database functionality with a simple BLAST command."""
        logger.debug("Testing BLAST database functionality")
        
        tmp_filename = None
        try:
            # Create temporary test query
            with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp_query:
                tmp_query.write(">test_seq\nACGTACGTACGTACGTACGT\n")
                tmp_query.flush()
                tmp_filename = tmp_query.name
            
            # Run test BLAST command
            cmd = [
                "blastn",
                "-task", "blastn-short",
                "-db", f'"{db_path}"',
                "-query", tmp_filename,
                "-outfmt", "6",
                "-max_target_seqs", "5"
            ]
            
            result = subprocess.run(cmd, text=True, capture_output=True, timeout=15)
            
            if result.returncode == 0:
                logger.debug("BLAST database test successful")
                return True
            else:
                # Check for non-fatal warnings
                error_msg = result.stderr.strip()
                if "memory map file error" in error_msg:
                    if re.search(r"Examining \d+ or more matches is recommended", error_msg):
                        logger.warning("Non-critical BLAST warning detected")
                        return True
                    
                    # Retry with blastdbcmd
                    return self._test_with_blastdbcmd(db_path)
                
                logger.error(f"BLAST test failed: {error_msg}")
                return False
                
        except subprocess.TimeoutExpired:
            logger.error("BLAST test timed out")
            return False
        except Exception as e:
            logger.error(f"Error testing database: {str(e)}")
            return False
        finally:
            if tmp_filename and os.path.exists(tmp_filename):
                try:
                    os.remove(tmp_filename)
                except OSError:
                    pass
    
    def _test_with_blastdbcmd(self, db_path: str) -> bool:
        """Fallback test using blastdbcmd."""
        try:
            cmd = ["blastdbcmd", "-db", f'"{db_path}"', "-info"]
            result = subprocess.run(cmd, text=True, capture_output=True, timeout=8)
            
            if result.returncode == 0:
                logger.debug("Database verified with blastdbcmd")
                return True
            else:
                logger.error("blastdbcmd verification failed")
                return False
                
        except Exception as e:
            logger.error(f"blastdbcmd test failed: {str(e)}")
            return False
    
    def _get_default_db_directory(self) -> str:
        """Get the default directory for storing BLAST databases."""
        if os.geteuid() == 0:  # Running as root
            return "/usr/local/share/ddprimer/blast_db"
        else:
            return os.path.join(Path.home(), ".ddprimer", "blast_db")
    
    def _get_organism_menu(self) -> str:
        """Generate the organism selection menu."""
        menu = "\nAvailable options:\n"
        menu += "0. Custom FASTA file\n"
        
        for i, (key, value) in enumerate(self.MODEL_ORGANISMS.items(), 1):
            menu += f"{i}. {value['name']} ({key})\n"
        
        menu += f"{len(self.MODEL_ORGANISMS) + 1}. Select from existing databases\n"
        menu += f"{len(self.MODEL_ORGANISMS) + 2}. Cancel\n"
        return menu
    
    def _create_from_model_organism(self) -> bool:
        """Create database from model organism selection."""
        logger.info("Creating database from model organism...")
        
        organism_key, organism_name, fasta_file = self.select_model_organism()
        
        if fasta_file is None:
            logger.info("Operation canceled.")
            return False
        
        if organism_key == 'existing_db':
            # User selected existing database
            return self.set_active_database(fasta_file)
        
        try:
            # Generate database name
            if organism_key is not None:
                organism_data = self.MODEL_ORGANISMS[organism_key]
                scientific_name = organism_data["name"].split(' (')[0]
                db_name = scientific_name.replace(' ', '_')
            else:
                db_name = None
            
            # Create the database
            db_path = self.create_database(fasta_file, db_name)
            
            # Set as active database
            if self.set_active_database(db_path):
                # Clean up genome file if from model organism
                if organism_key is not None:
                    self._cleanup_genome_file(fasta_file)
                return True
            else:
                return False
                
        except Exception as e:
            error_msg = f"Failed to create database from model organism: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return False
    
    def _create_from_custom_file(self) -> bool:
        """Create database from custom FASTA file."""
        logger.info("Creating database from custom FASTA file...")
        
        try:
            from ..utils import FileIO
            logger.info("\n>>> Please select a FASTA file for BLAST database creation <<<")
            fasta_file = FileIO.select_fasta_file("Select FASTA file for BLAST database creation")
            
            if fasta_file is None:
                logger.info("File selection canceled.")
                return False
            
            # Create database
            db_path = self.create_database(fasta_file)
            
            # Set as active database
            return self.set_active_database(db_path)
            
        except Exception as e:
            error_msg = f"Failed to create database from custom file: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return False
    
    def _select_existing_database(self) -> bool:
        """Select from existing databases."""
        logger.info("Selecting from existing databases...")
        
        databases = self.find_existing_databases()
        if not databases:
            logger.info("No existing databases found.")
            return False
        
        # Present selection menu
        logger.info("\nAvailable databases:")
        db_paths = list(databases.keys())
        
        for i, db_path in enumerate(db_paths, 1):
            display_name = databases[db_path]
            logger.info(f"{i}. {display_name}")
        
        logger.info(f"{len(databases) + 1}. Cancel")
        
        try:
            choice = input(f"\nEnter your choice [1-{len(databases) + 1}]: ")
            choice = int(choice)
            
            if 1 <= choice <= len(databases):
                selected_path = db_paths[choice - 1]
                return self.set_active_database(selected_path)
            else:
                logger.info("Selection canceled.")
                return False
                
        except (ValueError, KeyboardInterrupt):
            logger.info("Selection canceled.")
            return False
    
    def _handle_custom_file_selection(self) -> Tuple[Optional[str], Optional[str], Optional[str]]:
            """Handle custom file selection."""
            try:
                from ..utils import FileIO
                logger.info("\n>>> Please select a FASTA file <<<")
                fasta_file = FileIO.select_fasta_file("Select FASTA file")
                if fasta_file:
                    # Generate a database name from the filename
                    db_name = os.path.splitext(os.path.basename(fasta_file))[0]
                    return None, "Custom file", fasta_file
                else:
                    return None, None, None
            except Exception as e:
                logger.error(f"Error selecting custom file: {str(e)}")
                return None, None, None
    
    def _handle_existing_database_selection(self) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        """Handle existing database selection."""
        try:
            databases = self.find_existing_databases()
            if not databases:
                logger.info("No existing databases found.")
                return None, None, None
            
            # Present selection menu
            logger.info("\nAvailable databases:")
            db_paths = list(databases.keys())
            
            for i, db_path in enumerate(db_paths, 1):
                display_name = databases[db_path]
                logger.info(f"{i}. {display_name}")
            
            logger.info(f"{len(databases) + 1}. Cancel")
            
            choice = input(f"\nEnter your choice [1-{len(databases) + 1}]: ")
            choice = int(choice)
            
            if 1 <= choice <= len(databases):
                selected_path = db_paths[choice - 1]
                db_name = os.path.basename(selected_path).replace("_", " ")
                return 'existing_db', db_name, selected_path
            else:
                return None, None, None
                
        except (ValueError, KeyboardInterrupt):
            return None, None, None
        except Exception as e:
            logger.error(f"Error selecting existing database: {str(e)}")
            return None, None, None
    
    def _handle_model_organism_selection(self, choice: int) -> Tuple[Optional[str], Optional[str], Optional[str]]:
            """Handle model organism selection and download."""
            try:
                organism_key = list(self.MODEL_ORGANISMS.keys())[choice - 1]
                organism_data = self.MODEL_ORGANISMS[organism_key]
                organism_name = organism_data['name']
                
                # Create temporary directory for genome download
                temp_genome_dir = tempfile.mkdtemp(prefix="ddprimer_genomes_")
                logger.debug(f"Created temporary genome directory: {temp_genome_dir}")
                
                try:
                    # Fetch the genome
                    fasta_file = self._fetch_model_organism(organism_key, temp_genome_dir)
                    return organism_key, organism_name, fasta_file
                except Exception as e:
                    # Clean up temp directory on error
                    if os.path.exists(temp_genome_dir):
                        shutil.rmtree(temp_genome_dir)
                    raise e
                
            except Exception as e:
                logger.error(f"Error fetching model organism: {str(e)}")
                return None, None, None
    
    def _fetch_model_organism(self, organism_key: str, output_dir: str) -> str:
        """Fetch genome file for a model organism."""
        if organism_key not in self.MODEL_ORGANISMS:
            raise ValueError(f"Invalid organism key: {organism_key}")
        
        organism_data = self.MODEL_ORGANISMS[organism_key]
        url = organism_data['url']
        compressed = organism_data.get('compressed', False)
        
        # Determine output file path
        filename = os.path.basename(url)
        output_file = os.path.join(output_dir, filename)
        
        # Check if file already exists
        if compressed and filename.endswith('.gz'):
            uncompressed_file = output_file[:-3]
            if os.path.exists(uncompressed_file):
                logger.info(f"Genome file already exists: {uncompressed_file}")
                return uncompressed_file
        elif os.path.exists(output_file):
            logger.info(f"Genome file already exists: {output_file}")
            return output_file
        
        # Download the file
        logger.debug(f"\nDownloading {organism_data['name']} genome...")
        logger.info(f"Source: {url}")
        
        def progress_hook(blocknum, blocksize, totalsize):
            readsofar = blocknum * blocksize
            if totalsize > 0:
                percent = readsofar * 100 / totalsize
                sys.stdout.write(f"\rDownloading: {percent:.1f}% "
                               f"({readsofar/1024/1024:.1f} MB / {totalsize/1024/1024:.1f} MB)")
            else:
                sys.stdout.write(f"\rDownloading: {readsofar/1024/1024:.1f} MB")
            sys.stdout.flush()
        
        urllib.request.urlretrieve(url, output_file, progress_hook)
        sys.stdout.write("\n")
        
        # Extract if compressed
        if compressed and filename.endswith('.gz'):
            return self._extract_gzip(output_file)
        
        return output_file
    
    def _extract_gzip(self, gzip_file: str) -> str:
        """Extract a gzipped file."""
        output_file = gzip_file[:-3]  # Remove .gz extension
        logger.debug(f"Extracting {gzip_file} to {output_file}")
        
        try:
            with gzip.open(gzip_file, 'rb') as f_in:
                with open(output_file, 'wb') as f_out:
                    logger.info("Extracting compressed genome file...")
                    
                    chunk_size = 4 * 1024 * 1024  # 4 MB chunks
                    while True:
                        chunk = f_in.read(chunk_size)
                        if not chunk:
                            break
                        f_out.write(chunk)
            
            logger.debug("Extraction completed")
            return output_file
            
        except Exception as e:
            error_msg = f"Failed to extract {gzip_file}: {str(e)}"
            logger.error(error_msg)
            raise FileError(error_msg) from e
    
    def _cleanup_genome_file(self, file_path: str) -> None:
            """Clean up downloaded genome files and their temporary directory."""
            if not file_path:
                return
                
            try:
                # Get the directory containing the genome file
                genome_dir = os.path.dirname(file_path)
                
                # If it's a temporary directory (contains ddprimer_genomes_), remove the whole directory
                if "ddprimer_genomes_" in genome_dir:
                    logger.debug(f"Cleaning up temporary genome directory: {genome_dir}")
                    if os.path.exists(genome_dir):
                        shutil.rmtree(genome_dir)
                else:
                    # Fallback: just remove the specific file (for backwards compatibility)
                    if os.path.exists(file_path):
                        logger.debug(f"Cleaning up genome file: {file_path}")
                        os.remove(file_path)
                        
                        # Also remove .gz file if this was decompressed
                        if not file_path.endswith('.gz'):
                            gz_file = f"{file_path}.gz"
                            if os.path.exists(gz_file):
                                os.remove(gz_file)
                            
            except OSError as e:
                logger.warning(f"Error cleaning up genome files: {str(e)}")


# Convenience functions for backward compatibility and ease of use
def verify_blast_database() -> bool:
    """
    Verify that a valid BLAST database exists and can be accessed.
    
    Returns:
        True if database is valid, False otherwise
    """
    manager = BlastDatabaseManager()
    return manager.verify_database()


def setup_blast_database() -> bool:
    """
    Interactive setup for BLAST database when verification fails.
    
    Returns:
        True if database was successfully set up, False otherwise
    """
    manager = BlastDatabaseManager()
    return manager.setup_database_interactive()