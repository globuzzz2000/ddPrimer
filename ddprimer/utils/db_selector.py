#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Database Selector for ddPrimer pipeline.

Contains functionality for:
1. Finding existing BLAST databases in standard locations
2. Interactive database selection menus
3. Database path validation and verification
4. User-friendly database management

This module provides functionality for finding and selecting BLAST
databases from standard installation locations and user directories.
"""

import os
import logging
from pathlib import Path
from ..config import FileError

# Set up module logger
logger = logging.getLogger(__name__)


class DatabaseSelector:
    """
    Finds and manages BLAST database selection.
    
    This class provides methods for discovering existing BLAST databases
    in standard locations and presenting interactive selection menus
    for users to choose from available databases.
    
    Example:
        >>> databases = DatabaseSelector.find_existing_databases()
        >>> if databases:
        ...     selected = DatabaseSelector.select_database()
        ...     print(f"Selected: {selected}")
    """
    
    @staticmethod
    def find_existing_databases():
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
        logger.debug("=== DATABASE SEARCH DEBUG ===")
        databases = {}
        
        # Standard locations to search for databases
        possible_locations = [
            "/usr/local/share/ddprimer/blast_db",
            os.path.join(Path.home(), ".ddprimer", "blast_db")
        ]
        
        logger.debug(f"Searching in locations: {possible_locations}")
        
        for location in possible_locations:
            if os.path.exists(location):
                logger.debug(f"Searching for databases in {location}")
                
                try:
                    # Look for any .nhr files which indicate BLAST databases
                    db_files = [f for f in os.listdir(location) if f.endswith(".nhr")]
                    logger.debug(f"Found {len(db_files)} .nhr files in {location}")
                    
                    for db_file in db_files:
                        # Get database name without extension
                        db_name = os.path.splitext(db_file)[0]
                        db_path = os.path.join(location, db_name)
                        
                        # Use a more friendly display name
                        display_name = db_name.replace("_", " ")
                        
                        # Add to our dictionary
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
        
        logger.debug(f"Total databases found: {len(databases)}")
        logger.debug("=== END DATABASE SEARCH DEBUG ===")
        return databases
    
    @staticmethod
    def select_database():
        """
        Allow user to select from existing databases.
        
        Presents an interactive menu of available databases and allows
        the user to select one. Handles user input validation and
        provides clear feedback on selections.
        
        Returns:
            Selected database path or None if canceled/failed
            
        Raises:
            FileError: If database discovery fails
        """
        logger.debug("Starting database selection process")
        
        # Find existing databases
        try:
            databases = DatabaseSelector.find_existing_databases()
        except FileError:
            # Re-raise FileError without modification
            raise
        except Exception as e:
            error_msg = f"Error finding existing databases: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(error_msg) from e
        
        if not databases:
            logger.info("No existing databases found.")
            return None
            
        # Display menu of databases
        logger.info("\nAvailable BLAST databases:")
        db_paths = list(databases.keys())
        
        for i, db_path in enumerate(db_paths, 1):
            display_name = databases[db_path]
            logger.info(f"{i}. {display_name}")
            
        logger.info(f"{len(databases) + 1}. Cancel")
        
        try:
            # Get user selection
            choice = input(f"Enter your choice [1-{len(databases) + 1}]: ")
            
            try:
                choice = int(choice)
            except ValueError:
                logger.error("Invalid choice. Please enter a number.")
                return None
                
            if 1 <= choice <= len(databases):
                # Return the selected database path
                selected_path = db_paths[choice - 1]
                selected_name = databases[selected_path]
                logger.info(f"Selected database: {selected_name}")
                logger.debug(f"Selected database path: {selected_path}")
                return selected_path
            elif choice == len(databases) + 1:
                # User canceled
                logger.info("Database selection canceled.")
                return None
            else:
                logger.error("Invalid choice. Please enter a number within the range.")
                return None
                
        except KeyboardInterrupt:
            logger.info("\nOperation canceled by user.")
            return None
        except Exception as e:
            error_msg = f"Error during database selection: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return None