#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Model Organism Manager for ddPrimer pipeline.

Contains functionality for:
1. Model organism genome fetching and processing
2. Interactive organism selection menus
3. Genome file download with progress tracking
4. Compressed file extraction and management

This module provides functionality to fetch genome files for common
model organisms, supporting the BLAST database creation workflow
with automated genome retrieval from NCBI.
"""

import os
import sys
import logging
import urllib.request
import gzip
import time
from pathlib import Path
from ..config import FileError, ExternalToolError

# Set up module logger
logger = logging.getLogger(__name__)


class ModelOrganismManager:
    """
    Manages the fetching and processing of model organism genome files.
    
    This class provides methods for downloading genome files from NCBI
    for common model organisms, with support for compressed file handling
    and interactive organism selection.
    
    Attributes:
        MODEL_ORGANISMS: Dictionary of available model organisms with URLs
        
    Example:
        >>> key, name, path = ModelOrganismManager.select_model_organism()
        >>> if path:
        ...     print(f"Downloaded {name} genome to {path}")
    """
    
    # Dictionary of model organisms with their genome download URLs - ordered as requested
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
    
    @staticmethod
    def get_model_organism_menu():
        """
        Create a formatted menu of available model organisms.
        
        Generates a user-friendly menu string showing all available
        model organisms with their scientific names and common names.
        
        Returns:
            Formatted menu string for display to users
        """
        menu = "\nAvailable model organisms:\n"
        menu += "0. Custom FASTA file\n"
        
        # Add model organisms with index starting at 1
        for i, (key, value) in enumerate(ModelOrganismManager.MODEL_ORGANISMS.items(), 1):
            menu += f"{i}. {value['name']} ({key})\n"
            
        # Add option to select from existing databases
        menu += f"{len(ModelOrganismManager.MODEL_ORGANISMS) + 1}. Select from existing databases\n"
        menu += f"{len(ModelOrganismManager.MODEL_ORGANISMS) + 2}. Cancel\n"
        return menu
    
    @staticmethod
    def fetch_model_organism(key, output_dir=None):
        """
        Fetch and prepare the genome file for a model organism.
        
        Downloads the genome file for the specified organism from NCBI,
        handling compressed files and providing progress feedback.
        
        Args:
            key: Model organism key from MODEL_ORGANISMS dictionary
            output_dir: Directory to save the downloaded genome
            
        Returns:
            Path to the downloaded genome file
            
        Raises:
            ValueError: If the model organism key is invalid
            ExternalToolError: If the download fails
            FileError: If file operations fail
        """
        logger.debug("=== MODEL ORGANISM FETCH DEBUG ===")
        logger.debug(f"Fetching organism: {key}")
        
        # Check if model organism exists
        if key not in ModelOrganismManager.MODEL_ORGANISMS:
            error_msg = f"Invalid model organism key: {key}"
            logger.error(error_msg)
            raise ValueError(error_msg)
        
        organism_data = ModelOrganismManager.MODEL_ORGANISMS[key]
        organism_name = organism_data['name']
        url = organism_data['url']
        compressed = organism_data.get('compressed', False)
        
        logger.debug(f"Organism name: {organism_name}")
        logger.debug(f"Download URL: {url}")
        logger.debug(f"Compressed: {compressed}")
        
        # Create output directory if needed
        if output_dir is None:
            output_dir = os.path.join(Path.home(), "ddprimer_genomes")
        
        try:
            os.makedirs(output_dir, exist_ok=True)
        except OSError as e:
            error_msg = f"Failed to create output directory {output_dir}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(error_msg) from e
        
        # Determine the filename for the genome
        filename = os.path.basename(url)
        output_file = os.path.join(output_dir, filename)
        logger.debug(f"Output file: {output_file}")
        
        # Check if the file already exists
        if os.path.exists(output_file):
            logger.info(f"Genome file for {organism_name} already exists at {output_file}")
            
            # If compressed, extract it
            if compressed and filename.endswith('.gz'):
                uncompressed_file = output_file[:-3]  # Remove .gz extension
                if os.path.exists(uncompressed_file):
                    logger.info(f"Uncompressed genome file already exists at {uncompressed_file}")
                    logger.debug("=== END MODEL ORGANISM FETCH DEBUG ===")
                    return uncompressed_file
                else:
                    logger.info(f"Extracting compressed genome file {output_file}")
                    result = ModelOrganismManager._extract_gzip(output_file)
                    logger.debug("=== END MODEL ORGANISM FETCH DEBUG ===")
                    return result
            
            logger.debug("=== END MODEL ORGANISM FETCH DEBUG ===")
            return output_file
        
        # Download the genome file
        logger.info(f"Downloading genome for {organism_name}...")
        logger.info(f"Source URL: {url}")
        
        try:
            # Download with progress indicator
            def reporthook(blocknum, blocksize, totalsize):
                readsofar = blocknum * blocksize
                if totalsize > 0:
                    percent = readsofar * 100 / totalsize
                    sys.stdout.write(f"\rDownloading {filename}: {percent:.1f}% ({readsofar/1024/1024:.1f} MB / {totalsize/1024/1024:.1f} MB)")
                    sys.stdout.flush()
                else:
                    sys.stdout.write(f"\rDownloading {filename}: {readsofar/1024/1024:.1f} MB")
                    sys.stdout.flush()
            
            # Perform the download
            urllib.request.urlretrieve(url, output_file, reporthook)
            sys.stdout.write("\n")
            logger.debug(f"Download completed: {output_file}")
            
            # If compressed, extract it
            if compressed and filename.endswith('.gz'):
                logger.debug(f"Extracting compressed genome file {output_file}")
                result = ModelOrganismManager._extract_gzip(output_file)
                logger.debug("=== END MODEL ORGANISM FETCH DEBUG ===")
                return result
            
            logger.debug("=== END MODEL ORGANISM FETCH DEBUG ===")
            return output_file
            
        except Exception as e:
            error_msg = f"Failed to download genome file for {organism_name}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise ExternalToolError(error_msg, tool_name="wget/urllib") from e
    
    @staticmethod
    def _extract_gzip(gzip_file):
        """
        Extract a gzipped file.
        
        Decompresses a gzipped genome file with progress tracking
        and proper error handling for large files.
        
        Args:
            gzip_file: Path to the gzipped file
            
        Returns:
            Path to the extracted file
            
        Raises:
            FileError: If extraction fails or file is inaccessible
        """
        output_file = gzip_file[:-3]  # Remove .gz extension
        logger.debug(f"Extracting {gzip_file} to {output_file}")
        
        try:
            # Extract with progress indicator
            with gzip.open(gzip_file, 'rb') as f_in:
                # Get the file size
                f_in.seek(0, 2)  # Seek to the end
                file_size = f_in.tell()  # Get current position
                f_in.seek(0)  # Back to the beginning
                
                # Extract the file
                with open(output_file, 'wb') as f_out:
                    file_size_mb = file_size / (1024 * 1024)
                    logger.info(f"Extracting file (approximately {file_size_mb:.1f} MB)...")
                    
                    bytes_read = 0
                    chunk_size = 4 * 1024 * 1024  # 4 MB chunks
                    start_time = time.time()
                    
                    while True:
                        chunk = f_in.read(chunk_size)
                        if not chunk:
                            break
                        
                        f_out.write(chunk)
                        bytes_read += len(chunk)
                        
                        # Update progress every second
                        elapsed = time.time() - start_time
                        if elapsed > 1:
                            percent = bytes_read * 100 / file_size
                            sys.stdout.write(f"\rExtracting: {percent:.1f}% ({bytes_read/1024/1024:.1f} MB / {file_size_mb:.1f} MB)")
                            sys.stdout.flush()
                            start_time = time.time()
            
            sys.stdout.write("\n")
            logger.debug(f"Extraction completed: {output_file}")
            return output_file
            
        except (OSError, IOError, gzip.BadGzipFile) as e:
            error_msg = f"Failed to extract genome file {gzip_file}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(error_msg) from e
        except Exception as e:
            error_msg = f"Unexpected error extracting genome file {gzip_file}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(error_msg) from e
    
    @staticmethod
    def select_model_organism():
        """
        Present a menu to select a model organism or an existing database.
        
        Displays an interactive menu allowing users to choose from model
        organisms, custom files, or existing databases.
        
        Returns:
            Tuple of (organism_key, organism_name, file_path) or (None, None, None) if canceled.
            For existing database selection: ('existing_db', database_name, database_path)
            
        Raises:
            FileError: If file selection or database operations fail
        """
        logger.debug("Starting model organism selection")
        
        # Display the menu
        menu = ModelOrganismManager.get_model_organism_menu()
        logger.info(menu)
        
        try:
            # Get user selection
            choice = input(f"Enter your choice [0-{len(ModelOrganismManager.MODEL_ORGANISMS) + 2}]: ")
            
            try:
                choice = int(choice)
            except ValueError:
                logger.error("Invalid choice. Please enter a number.")
                return None, None, None
            
            # Process the selection
            if choice == 0:  # Custom file
                try:
                    from ..utils import FileIO
                    fasta_file = FileIO.select_fasta_file("Select FASTA file for BLAST database creation")
                    return None, "Custom file", fasta_file
                except Exception as e:
                    error_msg = f"Error selecting custom file: {str(e)}"
                    logger.error(error_msg)
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    return None, None, None
            elif choice == len(ModelOrganismManager.MODEL_ORGANISMS) + 1:  # Select from existing databases
                logger.info("Selecting from existing BLAST databases...")
                try:
                    from ..utils import DatabaseSelector
                    selected_db_path = DatabaseSelector.select_database()
                    
                    if selected_db_path is None:
                        logger.info("Database selection canceled.")
                        return None, None, None
                    
                    # Return special values to indicate an existing database was selected
                    db_name = os.path.basename(selected_db_path).replace("_", " ")
                    return 'existing_db', db_name, selected_db_path
                except Exception as e:
                    error_msg = f"Error selecting existing database: {str(e)}"
                    logger.error(error_msg)
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    return None, None, None
            elif choice == len(ModelOrganismManager.MODEL_ORGANISMS) + 2:  # Cancel
                return None, None, None
            elif 1 <= choice <= len(ModelOrganismManager.MODEL_ORGANISMS):
                # Get the organism key
                organism_key = list(ModelOrganismManager.MODEL_ORGANISMS.keys())[choice - 1]
                organism_data = ModelOrganismManager.MODEL_ORGANISMS[organism_key]
                organism_name = organism_data['name']
                
                # Create system directory for genome files
                if os.geteuid() == 0:  # Check if running as root
                    genome_dir = "/usr/local/share/ddprimer/genomes"
                else:
                    genome_dir = os.path.join(Path.home(), ".ddprimer", "genomes")
                
                try:
                    os.makedirs(genome_dir, exist_ok=True)
                except OSError as e:
                    error_msg = f"Failed to create genome directory {genome_dir}: {str(e)}"
                    logger.error(error_msg)
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    raise FileError(error_msg) from e
                
                # Fetch the genome file
                try:
                    file_path = ModelOrganismManager.fetch_model_organism(organism_key, genome_dir)
                    return organism_key, organism_name, file_path
                except (ValueError, ExternalToolError, FileError):
                    # Re-raise specific exceptions without modification
                    raise
                except Exception as e:
                    error_msg = f"Unexpected error fetching model organism: {str(e)}"
                    logger.error(error_msg)
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    raise FileError(error_msg) from e
            else:
                logger.error("Invalid choice. Please enter a number within the range.")
                return None, None, None
                
        except KeyboardInterrupt:
            logger.info("\nOperation canceled by user.")
            return None, None, None
        except (FileError, ValueError, ExternalToolError):
            # Re-raise specific exceptions without modification
            raise
        except Exception as e:
            error_msg = f"Error selecting model organism: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(error_msg) from e
        
    @staticmethod
    def cleanup_genome_file(file_path):
        """
        Delete a genome file after database creation.
        
        Removes downloaded genome files to free up disk space after
        successful database creation, including both compressed and
        uncompressed versions.
        
        Args:
            file_path: Path to the genome file to remove
        """
        if not file_path:
            logger.debug("No genome file path provided for cleanup")
            return
            
        try:
            if os.path.exists(file_path):
                logger.debug(f"Cleaning up genome file: {file_path}")
                os.remove(file_path)
                
                # If this was a decompressed file, also try to remove the original .gz file
                if not file_path.endswith('.gz'):
                    gz_file = f"{file_path}.gz"
                    if os.path.exists(gz_file):
                        logger.debug(f"Cleaning up compressed genome file: {gz_file}")
                        os.remove(gz_file)
                        
                logger.debug("Genome file cleanup completed")
            else:
                logger.debug(f"No genome file to clean up or file not found: {file_path}")
        except OSError as e:
            logger.warning(f"Error cleaning up genome file {file_path}: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
        except Exception as e:
            logger.warning(f"Unexpected error during genome file cleanup: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)