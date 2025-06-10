#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
General file selector module for ddPrimer pipeline.

Contains functionality for:
1. Finding existing BLAST databases in standard locations
2. Finding existing k-mer lists in standard locations
3. Interactive selection menus for files/databases
4. File path validation and verification
5. User-friendly file management

This module provides unified functionality for finding and selecting
various file types from standard installation locations and user directories.
"""

import os
import logging
from pathlib import Path
from ..config import FileError

# Set up module logger
logger = logging.getLogger(__name__)


class Selector:
    """
    General file and database selector for ddPrimer pipeline.
    
    This class provides methods for discovering existing files of various types
    in standard locations and presenting interactive selection menus for users.
    
    Example:
        >>> databases = Selector.find_blast_databases()
        >>> selected = Selector.select_from_menu(databases, "BLAST database")
        >>> kmer_lists = Selector.find_kmer_lists()
        >>> selected_kmer = Selector.select_from_menu(kmer_lists, "k-mer list")
    """
    
    @staticmethod
    def find_blast_databases():
        """
        Find all existing BLAST databases in standard locations.
        
        Searches common installation directories for BLAST database files
        (.nhr files) and builds a dictionary of available databases with
        user-friendly display names.
        
        Returns:
            Dictionary mapping database paths to display names
            
        Raises:
            FileError: If there are permission issues accessing directories
        """
        logger.debug("=== BLAST DATABASE SEARCH DEBUG ===")
        databases = {}
        
        # Standard locations to search for databases
        possible_locations = [
            "/usr/local/share/ddprimer/blast_db",
            os.path.join(Path.home(), ".ddprimer", "blast_db"),
            "/Library/Application Support/Blast_DBs"  # macOS location
        ]
        
        logger.debug(f"Searching for BLAST databases in: {possible_locations}")
        
        for location in possible_locations:
            if os.path.exists(location):
                logger.debug(f"Searching for databases in {location}")
                
                try:
                    # Look for .nhr files and subdirectories
                    for root, dirs, files in os.walk(location):
                        db_files = [f for f in files if f.endswith(".nhr")]
                        
                        for db_file in db_files:
                            # Get database name without extension
                            db_name = os.path.splitext(db_file)[0]
                            db_path = os.path.join(root, db_name)
                            
                            # Create friendly display name from path
                            rel_path = os.path.relpath(db_path, location)
                            display_name = rel_path.replace(os.sep, " / ").replace("_", " ")
                            
                            databases[db_path] = display_name
                            logger.debug(f"Found database: {display_name} at {db_path}")
                        
                except OSError as e:
                    error_msg = f"Permission error searching for databases in {location}: {str(e)}"
                    logger.error(error_msg)
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    raise FileError(error_msg) from e
                except Exception as e:
                    logger.warning(f"Error searching for databases in {location}: {str(e)}")
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
            else:
                logger.debug(f"Location does not exist: {location}")
        
        logger.debug(f"Total BLAST databases found: {len(databases)}")
        logger.debug("=== END BLAST DATABASE SEARCH DEBUG ===")
        return databases
    
    @staticmethod
    def find_kmer_lists():
        """
        Find all existing k-mer list files in standard locations, grouped by organism.
        
        Searches common directories for k-mer frequency list files
        (.list files) and groups them by organism name for easier selection.
        
        Returns:
            Dictionary mapping organism names to lists of k-mer file paths
            
        Raises:
            FileError: If there are permission issues accessing directories
        """
        logger.debug("=== K-MER LIST SEARCH DEBUG ===")
        kmer_lists_by_organism = {}
        
        # Standard locations to search for k-mer lists
        possible_locations = [
            os.path.join(Path.home(), ".ddprimer", "kmer_lists"),
            "/usr/local/share/ddprimer/kmer_lists",
            "./kmer_lists"  # Current directory
        ]
        
        logger.debug(f"Searching for k-mer lists in: {possible_locations}")
        
        for location in possible_locations:
            if os.path.exists(location):
                logger.debug(f"Searching for k-mer lists in {location}")
                
                try:
                    # Look for .list files
                    for root, dirs, files in os.walk(location):
                        list_files = [f for f in files if f.endswith(".list")]
                        
                        for list_file in list_files:
                            list_path = os.path.join(root, list_file)
                            
                            # Parse organism name from filename (e.g., "arabidopsis_thaliana_11.list")
                            filename_base = os.path.splitext(list_file)[0]
                            
                            # Split by underscore and extract organism name
                            parts = filename_base.split('_')
                            if len(parts) >= 3 and parts[-1].isdigit():
                                # Organism name is everything except the last part (k-mer size)
                                organism_name = '_'.join(parts[:-1]).replace('_', ' ')
                                kmer_size = parts[-1]
                                
                                # Group by organism
                                if organism_name not in kmer_lists_by_organism:
                                    kmer_lists_by_organism[organism_name] = []
                                
                                kmer_lists_by_organism[organism_name].append({
                                    'path': list_path,
                                    'size': kmer_size,
                                    'filename': list_file
                                })
                                
                                logger.debug(f"Found k-mer list: {organism_name} ({kmer_size}-mer) at {list_path}")
                        
                except OSError as e:
                    error_msg = f"Permission error searching for k-mer lists in {location}: {str(e)}"
                    logger.error(error_msg)
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    raise FileError(error_msg) from e
                except Exception as e:
                    logger.warning(f"Error searching for k-mer lists in {location}: {str(e)}")
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
            else:
                logger.debug(f"Location does not exist: {location}")
        
        # Sort k-mer files within each organism by size
        for organism in kmer_lists_by_organism:
            kmer_lists_by_organism[organism].sort(key=lambda x: int(x['size']))
        
        logger.debug(f"Total organisms with k-mer lists: {len(kmer_lists_by_organism)}")
        logger.debug("=== END K-MER LIST SEARCH DEBUG ===")
        return kmer_lists_by_organism
    
    @staticmethod
    def select_from_menu(items_dict, item_type="item"):
        """
        Present an interactive menu for item selection.
        
        Creates a numbered menu from the provided dictionary and allows
        the user to select an item. Handles user input validation and
        provides clear feedback on selections.
        
        Args:
            items_dict: Dictionary mapping item paths to display names
            item_type: Type of item for display purposes (e.g., "database", "k-mer list")
            
        Returns:
            Selected item path or None if canceled/failed
            
        Raises:
            FileError: If menu creation fails
        """
        logger.debug(f"Starting {item_type} selection process")
        
        if not items_dict:
            logger.info(f"No {item_type}s found.")
            return None
            
        # Display menu
        logger.info(f"\nAvailable {item_type}s:")
        item_paths = list(items_dict.keys())
        
        for i, item_path in enumerate(item_paths, 1):
            display_name = items_dict[item_path]
            logger.info(f"{i}. {display_name}")
            
        logger.info(f"{len(items_dict) + 1}. Cancel")
        
        try:
            # Get user selection
            choice = input(f"Enter your choice [1-{len(items_dict) + 1}]: ")
            
            try:
                choice = int(choice)
            except ValueError:
                logger.error("Invalid choice. Please enter a number.")
                return None
                
            if 1 <= choice <= len(items_dict):
                # Return the selected item path
                selected_path = item_paths[choice - 1]
                selected_name = items_dict[selected_path]
                logger.info(f"Selected {item_type}: {selected_name}")
                logger.debug(f"Selected {item_type} path: {selected_path}")
                return selected_path
            elif choice == len(items_dict) + 1:
                # User canceled
                logger.info(f"{item_type.capitalize()} selection canceled.")
                return None
            else:
                logger.error("Invalid choice. Please enter a number within the range.")
                return None
                
        except KeyboardInterrupt:
            logger.info("\nOperation canceled by user.")
            return None
        except Exception as e:
            error_msg = f"Error during {item_type} selection: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return None
    
    @staticmethod
    def select_blast_database():
        """
        Allow user to select from existing BLAST databases.
        
        Convenience method that combines database discovery and selection
        into a single operation.
        
        Returns:
            Selected database path or None if canceled/failed
            
        Raises:
            FileError: If database discovery fails
        """
        try:
            databases = Selector.find_blast_databases()
            return Selector.select_from_menu(databases, "BLAST database")
        except FileError:
            # Re-raise FileError without modification
            raise
        except Exception as e:
            error_msg = f"Error in BLAST database selection: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(error_msg) from e
    
    @staticmethod
    def select_kmer_list():
        """
        Allow user to select from existing k-mer lists grouped by organism.
        
        Convenience method that groups k-mer lists by organism and allows
        selection of all k-mer sizes for a chosen organism.
        
        Returns:
            List of selected k-mer list paths or None if canceled/failed
            
        Raises:
            FileError: If k-mer list discovery fails
        """
        try:
            kmer_lists_by_organism = Selector.find_kmer_lists()
            
            if not kmer_lists_by_organism:
                logger.info("No k-mer lists found.")
                return None
                
            # Display menu of organisms
            logger.info("\nAvailable k-mer lists by organism:")
            organisms = list(kmer_lists_by_organism.keys())
            
            for i, organism in enumerate(organisms, 1):
                sizes = [item['size'] for item in kmer_lists_by_organism[organism]]
                sizes_str = ', '.join(f"{size}-mer" for size in sizes)
                logger.info(f"{i}. {organism.title()} ({sizes_str})")
                
            logger.info(f"{len(organisms) + 1}. Cancel")
            
            try:
                # Get user selection
                choice = input(f"Enter your choice [1-{len(organisms) + 1}]: ")
                
                try:
                    choice = int(choice)
                except ValueError:
                    logger.error("Invalid choice. Please enter a number.")
                    return None
                    
                if 1 <= choice <= len(organisms):
                    # Return all k-mer files for the selected organism
                    selected_organism = organisms[choice - 1]
                    selected_files = [item['path'] for item in kmer_lists_by_organism[selected_organism]]
                    
                    sizes_str = ', '.join(f"{item['size']}-mer" for item in kmer_lists_by_organism[selected_organism])
                    logger.info(f"Selected {selected_organism.title()} k-mer lists: {sizes_str}")
                    
                    for path in selected_files:
                        logger.debug(f"Selected k-mer list path: {path}")
                    
                    return selected_files
                elif choice == len(organisms) + 1:
                    # User canceled
                    logger.info("K-mer list selection canceled.")
                    return None
                else:
                    logger.error("Invalid choice. Please enter a number within the range.")
                    return None
                    
            except KeyboardInterrupt:
                logger.info("\nOperation canceled by user.")
                return None
            except Exception as e:
                error_msg = f"Error during k-mer list selection: {str(e)}"
                logger.error(error_msg)
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                return None
                
        except FileError:
            # Re-raise FileError without modification
            raise
        except Exception as e:
            error_msg = f"Error in k-mer list selection: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(error_msg) from e