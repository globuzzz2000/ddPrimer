#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
General file selector module for ddPrimer pipeline.

Contains functionality for:
1. Finding existing BLAST databases in standard locations
2. Interactive selection menus for files/databases
3. File path validation and verification
4. User-friendly file management
5. Unified configuration management integration

This module provides unified functionality for finding and selecting
various file types from standard installation locations and user directories,
with integration to the unified configuration system.
"""

import os
import logging
from pathlib import Path
from typing import Optional
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
            item_type: Type of item for display purposes (e.g., "database")
            
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
        Display summary of available BLAST databases.
        
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
        
        # Show current selections
        try:
            db_config = Selector.load_database_config()
            if db_config:
                logger.info(f"\nCurrent BLAST Database: {Path(db_config).name}")
            else:
                logger.info("\nCurrent BLAST Database: None selected")
        except Exception as e:
            logger.debug(f"Error checking database config: {str(e)}")

    
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