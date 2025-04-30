#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Database Selector

A module for finding and selecting BLAST databases.
"""

import os
import logging
from pathlib import Path

class DatabaseSelector:
    """Finds and manages BLAST database selection."""
    
    @staticmethod
    def find_existing_databases(logger=None):
        """
        Find all existing BLAST databases in standard locations.
        
        Args:
            logger (logging.Logger, optional): Logger instance
            
        Returns:
            dict: Dictionary of {database_path: display_name}
        """
        if logger is None:
            logger = logging.getLogger("ddPrimer")
            
        databases = {}
        
        # Standard locations to search for databases
        possible_locations = [
            # System directory
            "/usr/local/share/ddprimer/blast_db",
            # User directory
            os.path.join(Path.home(), ".ddprimer", "blast_db")
        ]
        
        for location in possible_locations:
            if os.path.exists(location):
                logger.debug(f"Searching for databases in {location}")
                
                # Look for any .nhr files which indicate BLAST databases
                try:
                    db_files = [f for f in os.listdir(location) if f.endswith(".nhr")]
                    
                    for db_file in db_files:
                        # Get database name without extension
                        db_name = os.path.splitext(db_file)[0]
                        db_path = os.path.join(location, db_name)
                        
                        # Use a more friendly display name
                        display_name = db_name.replace("_", " ")
                        
                        # Add to our dictionary
                        databases[db_path] = display_name
                        logger.debug(f"Found database: {display_name} at {db_path}")
                        
                except Exception as e:
                    logger.warning(f"Error searching for databases in {location}: {str(e)}")
        
        return databases
    
    @staticmethod
    def select_database(logger=None):
        """
        Allow user to select from existing databases.
        
        Args:
            logger (logging.Logger, optional): Logger instance
            
        Returns:
            str: Selected database path or None if canceled
        """
        if logger is None:
            logger = logging.getLogger("ddPrimer")
            
        # Find existing databases
        databases = DatabaseSelector.find_existing_databases(logger)
        
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
                logger.info(f"Selected database: {databases[selected_path]}")
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
            logger.error(f"Error during database selection: {str(e)}")
            return None