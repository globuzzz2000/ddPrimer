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
6. Unified configuration management integration

This module provides unified functionality for finding and selecting
various file types from standard installation locations and user directories,
with integration to the unified configuration system.
"""

import os
import time
import logging
from pathlib import Path
from typing import Optional, List, Dict, Tuple
from ..config import FileError

# Set up module logger
logger = logging.getLogger(__name__)


class Selector:
    """
    General file and database selector for ddPrimer pipeline.
    
    This class provides methods for discovering existing files of various types
    in standard locations and presenting interactive selection menus for users.
    Now integrated with unified configuration system for persistent preferences.
    
    Example:
        >>> databases = Selector.find_blast_databases()
        >>> selected = Selector.select_from_menu(databases, "BLAST database")
        >>> kmer_lists = Selector.find_kmer_lists()
        >>> action, data = Selector.select_kmer_operation()
    """
    
    #############################################################################
    #                           GENERAL METHODS
    #############################################################################
    
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
                logger.debug(f"Selected {item_type}: {selected_name}")
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
    def show_available_resources():
        """
        Display summary of available BLAST databases and k-mer lists.
        
        Useful for debugging and user information about what resources
        are available on the system.
        """
        logger.info("\n=== Available ddPrimer Resources ===")
        
        # Show BLAST databases
        try:
            databases = Selector.find_blast_databases()
            if databases:
                logger.info(f"\nBLAST Databases ({len(databases)} found):")
                for db_path, display_name in databases.items():
                    logger.info(f"  - {display_name}")
                    logger.debug(f"    Path: {db_path}")
            else:
                logger.info("\nBLAST Databases: None found")
        except Exception as e:
            logger.warning(f"Error checking BLAST databases: {str(e)}")
        
        # Show k-mer lists
        try:
            kmer_lists = Selector.find_kmer_lists()
            if kmer_lists:
                logger.info(f"\nK-mer Lists ({len(kmer_lists)} prefixes):")
                for prefix, files in kmer_lists.items():
                    sizes = [item['size'] for item in files]
                    sizes_str = ', '.join(f"{size}-mer" for size in sizes)
                    logger.info(f"  - {prefix}: {sizes_str}")
                    for item in files:
                        logger.debug(f"    {item['filename']}: {item['path']}")
            else:
                logger.info("\nK-mer Lists: None found")
        except Exception as e:
            logger.warning(f"Error checking k-mer lists: {str(e)}")
        
        # Show current selections
        try:
            db_config = Selector.load_database_config()
            if db_config:
                logger.info(f"\nCurrent BLAST Database: {Path(db_config).name}")
            else:
                logger.info("\nCurrent BLAST Database: None selected")
        except Exception as e:
            logger.debug(f"Error checking database config: {str(e)}")
            
        try:
            kmer_info = Selector.get_kmer_selection_info()
            if kmer_info:
                logger.info(f"Current K-mer Selection: {kmer_info}")
            else:
                logger.info("Current K-mer Selection: None")
        except Exception as e:
            logger.debug(f"Error checking k-mer selection: {str(e)}")
        
        # Show tool availability
        glistmaker_available = Selector.check_glistmaker_availability()
        logger.info(f"\nGenomeTester4 (glistmaker): {'Available' if glistmaker_available else 'Not available'}")
        
        logger.info("=" * 40)
    
    #############################################################################
    #                           BLAST-SPECIFIC METHODS
    #############################################################################
    
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
    def select_blast_database():
        """
        Allow user to select from existing BLAST databases.
        
        Always presents selection menu rather than using saved preferences.
        
        Returns:
            Selected database path or None if canceled/failed
            
        Raises:
            FileError: If database discovery fails
        """
        try:
            # Always present selection menu (don't auto-load saved database)
            logger.info("\nPlease select a BLAST database:")
            databases = Selector.find_blast_databases()
            selected_db = Selector.select_from_menu(databases, "BLAST database")
            
            # Save the selection if user chose one
            if selected_db:
                Selector.save_database_config(selected_db)
                logger.debug(f"Selected BLAST database: {Path(selected_db).name}")
            
            return selected_db
            
        except FileError:
            # Re-raise FileError without modification
            raise
        except Exception as e:
            error_msg = f"Error in BLAST database selection: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(error_msg) from e
    
    @staticmethod
    def save_database_config(db_path: str):
        """
        Save database configuration using unified config system.
        
        Args:
            db_path: Path to the BLAST database
        """
        try:
            from ..config import Config
            Config.save_user_setting("blast_db_path", db_path)
            logger.debug(f"Saved database config to unified storage: {db_path}")
        except Exception as e:
            logger.warning(f"Failed to save database config: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)

    @staticmethod
    def load_database_config() -> Optional[str]:
        """
        Load database configuration from unified config system.
        
        Returns:
            Database path if found, None otherwise
        """
        try:
            from ..config import Config
            db_path = Config.load_user_setting("blast_db_path", None)
            
            if db_path and os.path.exists(db_path + ".nhr"):
                logger.debug(f"Loaded database config from unified storage: {db_path}")
                return db_path
            elif db_path:
                logger.debug(f"Saved database path no longer valid: {db_path}")
                return None
            else:
                logger.debug("No database config found in unified storage")
                return None
                
        except Exception as e:
            logger.warning(f"Error loading database config: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return None

    #############################################################################
    #                           K-MER SPECIFIC METHODS
    #############################################################################
    
    @staticmethod
    def find_kmer_lists():
        """
        Find all existing k-mer list files in standard locations, grouped by file prefix.
        
        Searches common directories for k-mer frequency list files
        (.list files) and groups them by file prefix for easier selection.
        
        Returns:
            Dictionary mapping file prefixes to lists of k-mer file info
            
        Raises:
            FileError: If there are permission issues accessing directories
        """
        logger.debug("=== K-MER LIST SEARCH DEBUG ===")
        kmer_lists_by_prefix = {}
        
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
                        logger.debug(f"Found .list files in {root}: {list_files}")
                        
                        for list_file in list_files:
                            list_path = os.path.join(root, list_file)
                            
                            # Parse prefix from filename (e.g., "thaliana_11.list" -> "thaliana")
                            filename_base = os.path.splitext(list_file)[0]
                            logger.debug(f"Processing file: {list_file}, base: {filename_base}")
                            
                            # Split by underscore and extract prefix
                            parts = filename_base.split('_')
                            logger.debug(f"Filename parts: {parts}")
                            
                            if len(parts) >= 2 and parts[-1].isdigit():
                                # Prefix is everything except the last part (k-mer size)
                                file_prefix = '_'.join(parts[:-1])
                                kmer_size = parts[-1]
                                
                                logger.debug(f"Parsed prefix: '{file_prefix}', k-mer size: {kmer_size}")
                                
                                # Group by prefix
                                if file_prefix not in kmer_lists_by_prefix:
                                    kmer_lists_by_prefix[file_prefix] = []
                                
                                kmer_lists_by_prefix[file_prefix].append({
                                    'path': list_path,
                                    'size': kmer_size,
                                    'filename': list_file,
                                    'folder': root
                                })
                                
                                logger.debug(f"Added k-mer list: {file_prefix} ({kmer_size}-mer) at {list_path}")
                            else:
                                logger.debug(f"Could not parse k-mer size from filename: {list_file}")
                        
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
        
        # Sort k-mer files within each prefix by size
        for prefix in kmer_lists_by_prefix:
            kmer_lists_by_prefix[prefix].sort(key=lambda x: int(x['size']))
        
        logger.debug(f"Total prefixes with k-mer lists: {len(kmer_lists_by_prefix)}")
        logger.debug(f"Found prefixes: {list(kmer_lists_by_prefix.keys())}")
        logger.debug("=== END K-MER LIST SEARCH DEBUG ===")
        return kmer_lists_by_prefix
    
    @staticmethod
    def check_glistmaker_availability() -> bool:
        """
        Check if glistmaker from GenomeTester4 is available.
        
        Returns:
            True if glistmaker is available, False otherwise
        """
        try:
            from ..utils import KmerGenerator
            generator = KmerGenerator()
            return generator.check_glistmaker_availability()
        except Exception as e:
            logger.debug(f"Error checking glistmaker availability: {str(e)}")
            return False

    @staticmethod
    def select_kmer_operation() -> Tuple[str, Optional[any]]:
        """
        Present menu for k-mer operations (generate new or select existing).
        
        Returns:
            Tuple of (action, data) where action is 'generate', 'select', or 'cancel'
            For 'generate': data is the FASTA file path
            For 'select': data is True (selection was successful)
            For 'cancel': data is None
        """
        try:
            # Get existing k-mer lists
            kmer_lists_by_prefix = Selector.find_kmer_lists()
            
            # Check glistmaker availability for generation
            glistmaker_available = Selector.check_glistmaker_availability()
            
            logger.info("\nK-mer List Options:")
            
            option_num = 1
            
            # Show existing k-mer lists if any
            if kmer_lists_by_prefix:
                logger.info(f"{option_num}. Select from existing k-mer lists")
                select_existing_option = option_num
                option_num += 1
                
                # Show what prefixes are available
                prefixes = list(kmer_lists_by_prefix.keys())
                logger.debug(f"   Available: {', '.join(prefixes)}")
            else:
                select_existing_option = None
                logger.info("   (No existing k-mer lists found)")
            
            # Show generation option if glistmaker is available
            if glistmaker_available:
                logger.info(f"{option_num}. Generate new k-mer lists from FASTA file")
                generate_option = option_num
                option_num += 1
            else:
                generate_option = None
                logger.info("(K-mer generation unavailable - glistmaker not found)")
                logger.info("Install GenomeTester4: conda install -c bioconda genometester4")
            
            logger.info(f"{option_num}. Cancel")
            cancel_option = option_num
            
            try:
                choice = input(f"Enter your choice [1-{option_num}]: ")
                choice = int(choice)
                
                # Handle selection of existing lists
                if select_existing_option and choice == select_existing_option:
                    selection_result = Selector.select_kmer_list()  # Now returns dict or None
                    if selection_result:
                        return 'select', selection_result
                    else:
                        return 'cancel', None
                
                # Handle generation of new lists
                elif generate_option and choice == generate_option:
                    from ..utils import FileIO
                    logger.info("\n>>> Please select a FASTA file for k-mer generation <<<")
                    fasta_file = FileIO.select_fasta_file("Select FASTA file for k-mer generation")
                    return 'generate', fasta_file
                
                # Handle cancel
                elif choice == cancel_option:
                    return 'cancel', None
                
                else:
                    logger.error("Invalid choice")
                    return 'cancel', None
                    
            except (ValueError, KeyboardInterrupt):
                logger.info("Operation canceled")
                return 'cancel', None
                
        except Exception as e:
            logger.error(f"Error in k-mer operation selection: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return 'cancel', None

    @staticmethod
    def select_kmer_list():
        """
        Allow user to select from existing k-mer lists grouped by file prefix.
        
        Returns:
            Dictionary with selection info if successful, None if canceled/failed
        """
        try:
            kmer_lists_by_prefix = Selector.find_kmer_lists()
            
            if not kmer_lists_by_prefix:
                logger.info("No k-mer lists found.")
                return None
                
            # Display menu of prefixes
            logger.info("\nAvailable k-mer lists by prefix:")
            prefixes = list(kmer_lists_by_prefix.keys())
            
            for i, prefix in enumerate(prefixes, 1):
                sizes = [item['size'] for item in kmer_lists_by_prefix[prefix]]
                sizes_str = ', '.join(f"{size}-mer" for size in sizes)
                logger.info(f"{i}. {prefix} ({sizes_str})")
                
            logger.info(f"{len(prefixes) + 1}. Cancel")
            
            try:
                # Get user selection
                choice = input(f"Enter your choice [1-{len(prefixes) + 1}]: ")
                
                try:
                    choice = int(choice)
                except ValueError:
                    logger.error("Invalid choice. Please enter a number.")
                    return None
                    
                if 1 <= choice <= len(prefixes):
                    # Get selected prefix info
                    selected_prefix = prefixes[choice - 1]
                    prefix_files = kmer_lists_by_prefix[selected_prefix]
                    
                    # Extract folder path from first file
                    first_file_path = prefix_files[0]['path']
                    folder_path = str(Path(first_file_path).parent)
                    
                    sizes_str = ', '.join(f"{item['size']}-mer" for item in prefix_files)
                    logger.debug(f"Selected {selected_prefix} k-mer lists: {sizes_str}")
                    logger.debug(f"K-mer folder: {folder_path}")
                    logger.debug(f"K-mer prefix: {selected_prefix}")
                    
                    # Save selection with folder path and prefix
                    Selector.save_kmer_selection(folder_path, selected_prefix)
                    
                    # Return selection info for backward compatibility
                    return {
                        'folder_path': folder_path,
                        'file_prefix': selected_prefix,
                        'files': [item['path'] for item in prefix_files],
                        'success': True
                    }
                elif choice == len(prefixes) + 1:
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
            error_msg = f"Error in k-mer list selection: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return None

    @staticmethod
    def save_kmer_selection(folder_path: str, file_prefix: str):
        """
        Save k-mer selection preference using folder path and file prefix.
        
        Args:
            folder_path: Path to folder containing k-mer files
            file_prefix: File prefix for k-mer files (e.g., "thaliana")
        """
        try:
            from ..config import Config
            kmer_config = {
                "folder_path": folder_path,
                "file_prefix": file_prefix,
                "selected_date": time.strftime('%Y-%m-%d %H:%M:%S')
            }
            Config.save_user_setting("preferred_kmer_config", kmer_config)
            logger.debug(f"Saved k-mer selection: folder={folder_path}, prefix={file_prefix}")
        except Exception as e:
            logger.warning(f"Failed to save k-mer selection: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)

    @staticmethod
    def load_kmer_selection() -> Optional[dict]:
        """
        Load k-mer selection preference with folder path and file prefix.
        
        Returns:
            Dictionary with folder_path and file_prefix, or None if not found
        """
        try:
            from ..config import Config
            kmer_config = Config.load_user_setting("preferred_kmer_config", None)
            
            if kmer_config and "folder_path" in kmer_config and "file_prefix" in kmer_config:
                # Validate that folder still exists
                folder_path = kmer_config["folder_path"]
                if os.path.exists(folder_path):
                    logger.debug(f"Loaded k-mer selection: folder={folder_path}, prefix={kmer_config['file_prefix']}")
                    return kmer_config
                else:
                    logger.debug(f"Saved k-mer folder no longer exists: {folder_path}")
                    return None
            else:
                logger.debug("No k-mer selection found in unified storage")
                return None
                
        except Exception as e:
            logger.warning(f"Error loading k-mer selection: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return None

    @staticmethod
    def clear_kmer_selection():
        """
        Clear saved k-mer selection from unified config system.
        
        Useful for resetting preferences or when k-mer files are no longer valid.
        """
        try:
            from ..config import Config
            Config.save_user_setting("preferred_kmer_config", None)
            logger.debug("Cleared k-mer selection from unified storage")
        except Exception as e:
            logger.warning(f"Failed to clear k-mer selection: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)

    @staticmethod
    def get_kmer_selection_info() -> Optional[str]:
        """
        Get human-readable information about current k-mer selection.
        
        Returns:
            String describing current k-mer selection, or None if none saved
        """
        try:
            kmer_config = Selector.load_kmer_selection()
            if kmer_config:
                folder_path = kmer_config.get("folder_path", "")
                file_prefix = kmer_config.get("file_prefix", "")
                date = kmer_config.get("selected_date", "unknown")
                
                return f"Prefix: {file_prefix} (folder: {Path(folder_path).name}) - selected on {date}"
                
            return None
            
        except Exception as e:
            logger.debug(f"Error getting k-mer selection info: {str(e)}")
            return None

    @staticmethod
    def validate_kmer_files(kmer_files: List[str]) -> List[str]:
        """
        Validate a list of k-mer files and return only valid ones.
        
        Args:
            kmer_files: List of k-mer file paths to validate
            
        Returns:
            List of valid k-mer file paths
        """
        valid_files = []
        
        for kmer_file in kmer_files:
            if os.path.exists(kmer_file) and os.path.getsize(kmer_file) > 0:
                # Basic validation - check if it's a .list file
                if kmer_file.endswith('.list'):
                    valid_files.append(kmer_file)
                    logger.debug(f"Validated k-mer file: {Path(kmer_file).name}")
                else:
                    logger.debug(f"Invalid k-mer file extension: {Path(kmer_file).name}")
            else:
                logger.debug(f"K-mer file missing or empty: {kmer_file}")
        
        return valid_files