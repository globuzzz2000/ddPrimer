#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Model Organism Manager

This module provides functionality to fetch genome files for common model organisms.
"""

import os
import sys
import logging
import urllib.request
import gzip
import time
from pathlib import Path

class ModelOrganismManager:
    """Manages the fetching and processing of model organism genome files."""
    
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
        
        Returns:
            str: Menu formatted as a string
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
    def fetch_model_organism(key, output_dir=None, logger=None):
        """
        Fetch and prepare the genome file for a model organism.
        
        Args:
            key (str): Model organism key (e.g., 'human', 'mouse')
            output_dir (str, optional): Directory to save the downloaded genome
            logger (logging.Logger, optional): Logger instance
            
        Returns:
            str: Path to the downloaded genome file
            
        Raises:
            ValueError: If the model organism key is invalid
            ConnectionError: If the download fails
        """
        # Use default logger if none provided
        if logger is None:
            logger = logging.getLogger("ddPrimer")
            
        # Check if model organism exists
        if key not in ModelOrganismManager.MODEL_ORGANISMS:
            raise ValueError(f"Invalid model organism key: {key}")
        
        organism_data = ModelOrganismManager.MODEL_ORGANISMS[key]
        organism_name = organism_data['name']
        url = organism_data['url']
        compressed = organism_data.get('compressed', False)
        
        # Create output directory if needed
        if output_dir is None:
            # Use a standard location in the user's home directory
            output_dir = os.path.join(Path.home(), "ddprimer_genomes")
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Determine the filename for the genome
        filename = os.path.basename(url)
        output_file = os.path.join(output_dir, filename)
        
        # Check if the file already exists
        if os.path.exists(output_file):
            logger.info(f"Genome file for {organism_name} already exists at {output_file}")
            
            # If compressed, extract it
            if compressed and filename.endswith('.gz'):
                uncompressed_file = output_file[:-3]  # Remove .gz extension
                if os.path.exists(uncompressed_file):
                    logger.info(f"Uncompressed genome file already exists at {uncompressed_file}")
                    return uncompressed_file
                else:
                    logger.info(f"Extracting compressed genome file {output_file}")
                    return ModelOrganismManager._extract_gzip(output_file, logger)
            
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
                return ModelOrganismManager._extract_gzip(output_file, logger)
            
            return output_file
            
        except Exception as e:
            logger.error(f"Error downloading genome file: {str(e)}")
            raise ConnectionError(f"Failed to download genome file: {str(e)}")
    
    @staticmethod
    def _extract_gzip(gzip_file, logger=None):
        """
        Extract a gzipped file.
        
        Args:
            gzip_file (str): Path to the gzipped file
            logger (logging.Logger, optional): Logger instance
            
        Returns:
            str: Path to the extracted file
            
        Raises:
            IOError: If extraction fails
        """
        # Use default logger if none provided
        if logger is None:
            logger = logging.getLogger("ddPrimer")
            
        output_file = gzip_file[:-3]  # Remove .gz extension
        
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
            
        except Exception as e:
            logger.error(f"Error extracting genome file: {str(e)}")
            raise IOError(f"Failed to extract genome file: {str(e)}")
    
    @staticmethod
    def select_model_organism(logger=None):
        """
        Present a menu to select a model organism or an existing database.
        
        Args:
            logger (logging.Logger, optional): Logger instance
            
        Returns:
            tuple: (organism_key, organism_name, file_path) or (None, None, None) if canceled
                  For existing database selection: ('existing_db', database_name, None)
        """
        # Use default logger if none provided
        if logger is None:
            logger = logging.getLogger("ddPrimer")
            
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
                from ..utils import FileIO
                try:
                    fasta_file = FileIO.select_fasta_file("Select FASTA file for BLAST database creation")
                    return None, "Custom file", fasta_file
                except Exception as e:
                    logger.error(f"Error selecting custom file: {str(e)}")
                    return None, None, None
            elif choice == len(ModelOrganismManager.MODEL_ORGANISMS) + 1:  # Select from existing databases
                logger.info("Selecting from existing BLAST databases...")
                from ..utils import DatabaseSelector
                selected_db_path = DatabaseSelector.select_database(logger)
                
                if selected_db_path is None:
                    logger.info("Database selection canceled.")
                    return None, None, None
                
                # Return special values to indicate an existing database was selected
                # We set file_path to None since we're not downloading a genome file
                db_name = os.path.basename(selected_db_path).replace("_", " ")
                return 'existing_db', db_name, selected_db_path
            elif choice == len(ModelOrganismManager.MODEL_ORGANISMS) + 2:  # Cancel
                return None, None, None
            elif 1 <= choice <= len(ModelOrganismManager.MODEL_ORGANISMS):
                # Get the organism key
                organism_key = list(ModelOrganismManager.MODEL_ORGANISMS.keys())[choice - 1]
                organism_data = ModelOrganismManager.MODEL_ORGANISMS[organism_key]
                organism_name = organism_data['name']
                
                # Create system directory for genome files - use /usr/local/share if root, otherwise user home
                if os.geteuid() == 0:  # Check if running as root
                    genome_dir = "/usr/local/share/ddprimer/genomes"
                else:
                    genome_dir = os.path.join(Path.home(), ".ddprimer", "genomes")
                
                os.makedirs(genome_dir, exist_ok=True)
                
                # Fetch the genome file
                file_path = ModelOrganismManager.fetch_model_organism(organism_key, genome_dir, logger)
                
                return organism_key, organism_name, file_path
            else:
                logger.error("Invalid choice. Please enter a number within the range.")
                return None, None, None
                
        except KeyboardInterrupt:
            logger.info("\nOperation canceled by user.")
            return None, None, None
        except Exception as e:
            logger.error(f"Error selecting model organism: {str(e)}")
            return None, None, None
        
    @staticmethod
    def cleanup_genome_file(file_path, logger=None):
        """
        Delete a genome file after database creation.
        
        Args:
            file_path (str): Path to the genome file
            logger (logging.Logger, optional): Logger instance
        """
        if logger is None:
            logger = logging.getLogger("ddPrimer")
        
        try:
            if file_path and os.path.exists(file_path):
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
        except Exception as e:
            logger.warning(f"Error cleaning up genome file: {str(e)}")