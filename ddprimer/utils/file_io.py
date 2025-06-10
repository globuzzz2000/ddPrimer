#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File I/O module for ddPrimer pipeline.

Contains functionality for:
1. Cross-platform file selection with GUI and CLI support
2. FASTA file reading and writing operations
3. Excel file formatting with comprehensive styling
4. Sequence table loading from CSV and Excel formats
5. Temporary directory management and cleanup

This module provides consistent interfaces for file operations including
file selection, reading/writing, and format conversions across different
platforms and environments.
"""

import os
import sys
import platform
import logging
import pandas as pd
import contextlib
import tempfile
import shutil
from ..config import FileSelectionError, FileFormatError, FileError

# Set up module logger
logger = logging.getLogger(__name__)

# Optional import for Excel formatting
try:
    import openpyxl
    from openpyxl.styles import PatternFill, Alignment, Font, Border, Side
    from openpyxl.utils import get_column_letter
    HAS_OPENPYXL = True
except ImportError:
    HAS_OPENPYXL = False

# Optional import for wxPython file dialogs
try:
    import wx
    HAS_WX = True
except ImportError:
    HAS_WX = False


@contextlib.contextmanager
def _silence_cocoa_stderr():
    """
    Temporarily redirect C-level stderr to /dev/null (macOS IMK chatter).
    
    This context manager suppresses macOS Input Method Kit warning messages
    that can clutter the console output during GUI operations.
    
    Yields:
        None - Context manager for silencing stderr on macOS
    """
    if platform.system() != "Darwin":
        yield
        return
        
    import os, sys
    old_fd = os.dup(2)
    devnull = os.open(os.devnull, os.O_WRONLY)
    os.dup2(devnull, 2)
    os.close(devnull)
    try:
        yield
    finally:
        os.dup2(old_fd, 2)
        os.close(old_fd)


class FileIO:
    """
    Handles all file I/O operations with consistent interfaces.
    
    This class provides cross-platform file operations including GUI and CLI
    file selection, FASTA file processing, Excel formatting, and sequence
    table loading. It automatically detects the environment and chooses
    appropriate methods for file operations.
    
    Attributes:
        use_cli: Whether to use CLI mode instead of GUI
        is_macos: Whether running on macOS
        has_pyobjc: Whether PyObjC is available on macOS
        
    Example:
        >>> sequences = FileIO.load_fasta("genome.fasta")
        >>> file_path = FileIO.select_fasta_file("Select genome file")
        >>> FileIO.save_results(df, "output_dir", "input.fasta")
    """
    
    # Default to GUI mode unless explicitly detected as headless
    use_cli = False

    # Check for macOS and PyObjC availability
    is_macos = platform.system() == "Darwin"
    has_pyobjc = False
    
    if is_macos:
        try:
            import Foundation
            import AppKit
            has_pyobjc = True
        except ImportError:
            has_pyobjc = False

    # Detect headless environments explicitly with refined logic
    if platform.system() == "Linux":
        if not os.environ.get("DISPLAY", ""):
            use_cli = True
            reason = "No DISPLAY environment variable on Linux."
        else:
            reason = "DISPLAY variable found on Linux."
    elif platform.system() == "Windows":
        if not sys.stdout.isatty():
            use_cli = True
            reason = "Non-interactive terminal on Windows."
        else:
            reason = "Interactive terminal on Windows."
    elif platform.system() == "Darwin":
        # Force GUI mode by default on macOS unless explicitly headless
        use_cli = False
        reason = "Defaulting to GUI mode on macOS."
    else:
        reason = "Unknown platform, defaulting to GUI mode."

    # Initialize last directory with user's home directory
    _last_directory = None
    
    # Shared wxPython app instance to manage file dialogs
    _wx_app = None
    
    @classmethod
    def initialize_wx_app(cls):
        """
        Initialize the wxPython app if it doesn't exist yet.
        
        Creates a wxPython application instance for file dialogs,
        with proper error handling and logging.
        
        Raises:
            FileSelectionError: If wxPython initialization fails
        """
        if HAS_WX and cls._wx_app is None:
            try:
                with _silence_cocoa_stderr():
                    cls._wx_app = wx.App(False)
                    logger.debug("wxPython app initialized successfully")
            except Exception as e:
                error_msg = f"Failed to initialize wxPython app: {str(e)}"
                logger.error(error_msg)
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                raise FileSelectionError(error_msg) from e
    
    @classmethod
    def hide_app(cls):
        """
        Hide the wxPython app without destroying it.
        
        This is safer on macOS and avoids segmentation faults by hiding
        the application from the dock rather than destroying it.
        """
        if cls.is_macos and cls.has_pyobjc and cls._wx_app is not None:
            try:
                import AppKit
                # Hide from Dock
                NSApplication = AppKit.NSApplication.sharedApplication()
                NSApplication.setActivationPolicy_(AppKit.NSApplicationActivationPolicyAccessory)
                logger.debug("App hidden from macOS dock")
            except Exception as e:
                logger.warning(f"Error hiding app from dock: {str(e)}")
                logger.debug(f"Error details: {str(e)}", exc_info=True)

    @classmethod
    def mark_selection_complete(cls):
        """
        Mark that all file selections are complete and hide the app from the dock.
        
        Should be called when file selection operations are finished to clean
        up the GUI application properly.
        """
        logger.debug("File selection process marked as complete")
        cls.hide_app()

    @classmethod
    def get_last_directory(cls):
        """
        Get the last directory used, loading from persistent storage if needed.
        
        Loads the last used directory from configuration file to provide
        a better user experience by starting file dialogs in the most
        recently used location.
        
        Returns:
            Path to the last directory used, or home directory as fallback
        """
        if cls._last_directory is None:
            logger.debug("=== LAST DIRECTORY LOADING DEBUG ===")
            config_dir = os.path.join(os.path.expanduser("~"), ".ddPrimer")
            
            try:
                os.makedirs(config_dir, exist_ok=True)
            except OSError as e:
                logger.warning(f"Failed to create config directory {config_dir}: {str(e)}")
                cls._last_directory = os.path.expanduser("~")
                logger.debug("=== END LAST DIRECTORY LOADING DEBUG ===")
                return cls._last_directory
                
            config_file = os.path.join(config_dir, "last_directory.txt")
            
            if os.path.exists(config_file):
                try:
                    with open(config_file, 'r') as f:
                        last_dir = f.read().strip()
                    if os.path.isdir(last_dir):
                        cls._last_directory = last_dir
                        logger.debug(f"Loaded last directory from config: {last_dir}")
                    else:
                        cls._last_directory = os.path.expanduser("~")
                        logger.debug(f"Last directory no longer exists, using home: {cls._last_directory}")
                except (OSError, IOError) as e:
                    logger.warning(f"Error reading last directory config: {str(e)}")
                    cls._last_directory = os.path.expanduser("~")
                except Exception as e:
                    logger.warning(f"Unexpected error reading last directory config: {str(e)}")
                    logger.debug(f"Error details: {str(e)}", exc_info=True)
                    cls._last_directory = os.path.expanduser("~")
            else:
                cls._last_directory = os.path.expanduser("~")
                logger.debug(f"No last directory config found, using home: {cls._last_directory}")
            
            logger.debug("=== END LAST DIRECTORY LOADING DEBUG ===")
        
        return cls._last_directory
    
    @classmethod
    def save_last_directory(cls, directory):
        """
        Save the last directory to persistent storage.
        
        Persists the directory path to configuration file for future use,
        improving user experience in subsequent file selections.
        
        Args:
            directory: Path to save as last used directory
        """
        if not directory or not isinstance(directory, str):
            logger.debug("Invalid directory provided to save_last_directory")
            return
            
        if not os.path.isdir(directory):
            logger.debug(f"Directory does not exist, not saving: {directory}")
            return
            
        cls._last_directory = directory
        
        # Save to config file
        config_dir = os.path.join(os.path.expanduser("~"), ".ddPrimer")
        
        try:
            os.makedirs(config_dir, exist_ok=True)
            config_file = os.path.join(config_dir, "last_directory.txt")
            
            with open(config_file, 'w') as f:
                f.write(directory)
            logger.debug(f"Saved last directory to config: {directory}")
        except (OSError, IOError) as e:
            logger.warning(f"Error saving last directory config: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
        except Exception as e:
            logger.warning(f"Unexpected error saving last directory config: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)

    @classmethod
    def normalize_filetypes(cls, filetypes):
        """
        Normalize file types for different platforms.
        
        Converts file type specifications to platform-appropriate formats,
        handling differences between macOS, Windows, and Linux file dialogs.
        
        Args:
            filetypes: List of (description, extension) tuples
            
        Returns:
            Normalized file types list for the current platform
        """
        logger.debug(f"Normalizing {len(filetypes) if filetypes else 0} file types for platform {platform.system()}")
        
        # Define "All Files" for different platforms
        if cls.is_macos or platform.system() == "Linux":
            all_files = ("All Files", "*")
        else:  # Windows
            all_files = ("All Files", "*.*")
            
        # Start with all files as a default
        normalized = [all_files]
        
        # Add specific file types
        if filetypes:
            for desc, ext in filetypes:
                if ext == "*" or ext == "*.*":
                    continue  # Skip generic "all files" entries
                
                if cls.is_macos:
                    # macOS: extension without wildcards (e.g., "txt" not "*.txt")
                    clean_ext = ext.lstrip("*.")
                    
                    # Handle compound extensions for macOS
                    if "." in clean_ext:
                        # For compound extensions like vcf.gz, use just the last part
                        base_ext = clean_ext.split(".")[-1]
                        if base_ext:
                            normalized.insert(0, (desc, base_ext))
                    elif clean_ext:
                        normalized.insert(0, (desc, clean_ext))
                else:
                    # Windows/Linux: ensure wildcard format (e.g., "*.txt")
                    if not ext.startswith("*.") and ext != "*":
                        clean_ext = f"*.{ext.lstrip('.')}"
                    else:
                        clean_ext = ext
                    normalized.insert(0, (desc, clean_ext))
        
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Normalized file types:")
            for desc, ext in normalized:
                logger.debug(f"  {desc}: {ext}")
                    
        return normalized
        
    @classmethod
    def select_file(cls, prompt, filetypes):
        """
        Show a file dialog or prompt for a file path if CLI mode is enabled.
        
        Provides cross-platform file selection with automatic fallback from
        GUI to CLI mode when necessary. Handles file validation and user
        experience improvements like remembering last directory.
        
        Args:
            prompt: Text to display in the file dialog
            filetypes: List of file type tuples for the dialog
            
        Returns:
            Selected file path
            
        Raises:
            FileSelectionError: If file selection fails or user cancels
        """
        logger.debug("=== FILE SELECTION DEBUG ===")
        
        # Get the last directory from our persistent storage
        last_directory = cls.get_last_directory()
        logger.debug(f"Starting file selection from directory: {last_directory}")

        if cls.use_cli:
            logger.debug(f"{prompt} (CLI mode)")
            try:
                file_path = input(f"{prompt}: ").strip()
                if not file_path:
                    error_msg = "No file path provided in CLI mode"
                    logger.error(error_msg)
                    raise FileSelectionError(error_msg)
                    
                # Validate the path exists
                if not os.path.exists(file_path):
                    error_msg = f"File not found: {file_path}"
                    logger.error(error_msg)
                    raise FileSelectionError(error_msg)
                
                logger.debug(f"CLI file selection successful: {file_path}")
                logger.debug("=== END FILE SELECTION DEBUG ===")
                return file_path
            except KeyboardInterrupt:
                error_msg = "File selection canceled by user"
                logger.error(error_msg)
                raise FileSelectionError(error_msg)
            except Exception as e:
                error_msg = f"CLI file selection failed: {str(e)}"
                logger.error(error_msg)
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                raise FileSelectionError(error_msg) from e

        # wxPython implementation
        if not cls.use_cli and HAS_WX:
            try:
                # Initialize the wx.App if it doesn't exist yet
                cls.initialize_wx_app()
                
                valid_filetypes = cls.normalize_filetypes(filetypes)
                # Exclude the generic 'All Files' filter so users see only specific types
                specific_types = [ft for ft in valid_filetypes if not ft[0].startswith("All Files")]
                wildcard = "|".join([f"{desc} ({ext})|{ext}" for desc, ext in specific_types]) or "All Files (*.*)|*.*"

                with _silence_cocoa_stderr():
                    dlg = wx.FileDialog(
                        None,
                        prompt,
                        wildcard=wildcard,
                        style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
                    )
                    dlg.SetDirectory(last_directory)
                    # Ensure the first specific filter is pre‑selected
                    dlg.SetFilterIndex(0)

                    if dlg.ShowModal() == wx.ID_OK:
                        file_path = dlg.GetPath()
                    else:
                        file_path = ""
                    dlg.Destroy()

                if not file_path:
                    error_msg = "No file was selected in the dialog"
                    logger.error(error_msg)
                    raise FileSelectionError(error_msg)

                # Update the last directory
                cls.save_last_directory(os.path.dirname(file_path))
                logger.debug(f"GUI file selection successful: {file_path}")
                logger.debug("=== END FILE SELECTION DEBUG ===")
                return file_path

            except FileSelectionError:
                # Re-raise FileSelectionError without modification
                raise
            except Exception as e:
                logger.warning(f"wxPython file selection failed: {str(e)}")
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                # Fall back to CLI mode if wxPython fails
                cls.use_cli = True
        
        # If wxPython isn't available or also failed, use CLI
        logger.warning("Falling back to CLI mode for file selection")
        return cls.select_file(prompt, filetypes)
        
    @classmethod
    def select_files(cls, prompt, filetypes):
        """
        Show a file dialog for multiple files or prompt if CLI mode is enabled.
        
        Provides multi-file selection with the same cross-platform support
        and fallback mechanisms as single file selection.
        
        Args:
            prompt: Text to display in the file dialog
            filetypes: List of file type tuples for the dialog
            
        Returns:
            List of selected file paths
            
        Raises:
            FileSelectionError: If file selection fails or user cancels
        """
        logger.debug("=== MULTIPLE FILE SELECTION DEBUG ===")
        
        # Get the last directory from our persistent storage
        last_directory = cls.get_last_directory()
        logger.debug(f"Starting multiple file selection from directory: {last_directory}")

        if cls.use_cli:
            logger.info(f"{prompt} (CLI mode)")
            try:
                paths = input(f"{prompt} (paths separated by spaces): ").strip()
                file_paths = paths.split() if paths else []
                if not file_paths:
                    error_msg = "No file paths provided in CLI mode"
                    logger.error(error_msg)
                    raise FileSelectionError(error_msg)
                    
                # Validate the paths exist
                for path in file_paths:
                    if not os.path.exists(path):
                        error_msg = f"File not found: {path}"
                        logger.error(error_msg)
                        raise FileSelectionError(error_msg)
                
                logger.debug(f"CLI multiple file selection successful: {len(file_paths)} files")
                logger.debug("=== END MULTIPLE FILE SELECTION DEBUG ===")
                return file_paths
            except KeyboardInterrupt:
                error_msg = "Multiple file selection canceled by user"
                logger.error(error_msg)
                raise FileSelectionError(error_msg)
            except Exception as e:
                error_msg = f"CLI multiple file selection failed: {str(e)}"
                logger.error(error_msg)
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                raise FileSelectionError(error_msg) from e

        # wxPython implementation
        if not cls.use_cli and HAS_WX:
            try:
                # Initialize the wx.App if it doesn't exist yet
                cls.initialize_wx_app()
                
                valid_filetypes = cls.normalize_filetypes(filetypes)
                # Exclude the generic 'All Files' filter so users see only specific types
                specific_types = [ft for ft in valid_filetypes if not ft[0].startswith("All Files")]
                wildcard = "|".join([f"{desc} ({ext})|{ext}" for desc, ext in specific_types]) or "All Files (*.*)|*.*"

                with _silence_cocoa_stderr():
                    dlg = wx.FileDialog(
                        None,
                        prompt,
                        wildcard=wildcard,
                        style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_MULTIPLE
                    )
                    dlg.SetDirectory(last_directory)
                    # Ensure the first specific filter is pre‑selected
                    dlg.SetFilterIndex(0)

                    file_paths = []
                    if dlg.ShowModal() == wx.ID_OK:
                        file_paths = dlg.GetPaths()
                    dlg.Destroy()

                if not file_paths:
                    error_msg = "No files were selected in the dialog"
                    logger.error(error_msg)
                    raise FileSelectionError(error_msg)

                # Update the last directory
                cls.save_last_directory(os.path.dirname(file_paths[0]))
                logger.debug(f"GUI multiple file selection successful: {len(file_paths)} files")
                logger.debug("=== END MULTIPLE FILE SELECTION DEBUG ===")
                return file_paths

            except FileSelectionError:
                # Re-raise FileSelectionError without modification
                raise
            except Exception as e:
                logger.warning(f"wxPython multiple file selection failed: {str(e)}")
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                # Fall back to CLI mode if wxPython fails
                cls.use_cli = True

        # If wxPython isn't available or also failed, use CLI
        logger.warning("Falling back to CLI mode for multiple file selection")
        return cls.select_files(prompt, filetypes)
        
    @staticmethod
    def load_fasta(filepath):
        """
        Load sequences from a FASTA file into a dictionary.
        
        Efficiently parses FASTA files with memory optimization for large
        genomes. Handles validation and error reporting comprehensively.
        
        Args:
            filepath: Path to the FASTA file
            
        Returns:
            Dictionary mapping sequence headers to sequences
            
        Raises:
            FileError: If the FASTA file doesn't exist
            FileFormatError: If there's an error parsing the FASTA file
        """
        logger.debug("=== FASTA LOADING DEBUG ===")
        sequences = {}
        name = None
        seq_chunks = []
        
        logger.debug(f"Loading FASTA file: {filepath}")
        
        # Validate file exists
        if not os.path.exists(filepath):
            error_msg = f"FASTA file not found: {filepath}"
            logger.error(error_msg)
            raise FileError(error_msg)
            
        try:
            with open(filepath, 'r') as f:
                line_count = 0
                for line in f:
                    line_count += 1
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if name:
                            sequences[name] = "".join(seq_chunks).upper()
                            logger.debug(f"Loaded sequence: {name} ({len(sequences[name]):,} bp)")
                        name = line[1:].split()[0]
                        seq_chunks = []
                        logger.debug(f"Started new sequence: {name}")
                    else:
                        seq_chunks.append(line)
                        
                # Handle the last sequence
                if name:
                    sequences[name] = "".join(seq_chunks).upper()
                    logger.debug(f"Loaded final sequence: {name} ({len(sequences[name]):,} bp)")
                    
        except (OSError, IOError) as e:
            error_msg = f"Error reading FASTA file {os.path.abspath(filepath)}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END FASTA LOADING DEBUG ===")
            raise FileError(error_msg) from e
        except Exception as e:
            error_msg = f"Error parsing FASTA file {os.path.abspath(filepath)}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END FASTA LOADING DEBUG ===")
            raise FileFormatError(error_msg) from e
        
        logger.debug(f"Successfully loaded {len(sequences)} sequences from FASTA file")
        logger.debug("=== END FASTA LOADING DEBUG ===")
        return sequences
    
    @staticmethod
    def save_fasta(sequences, filepath):
        """
        Save sequences to a FASTA file.
        
        Writes sequences in proper FASTA format with 60-character line wrapping
        for optimal readability and compatibility.
        
        Args:
            sequences: Dictionary of sequence ID to sequence
            filepath: Path to save the FASTA file
            
        Raises:
            FileFormatError: If there's an error writing the FASTA file
        """
        logger.info(f"Saving {len(sequences)} sequences to FASTA file: {filepath}")
        
        if not sequences:
            logger.warning("No sequences provided to save_fasta")
            
        try:
            with open(filepath, 'w') as f:
                for seq_id, sequence in sequences.items():
                    if not seq_id or not sequence:
                        logger.warning(f"Skipping invalid sequence: id='{seq_id}', seq_len={len(sequence) if sequence else 0}")
                        continue
                        
                    f.write(f">{seq_id}\n")
                    # Write sequence in chunks of 60 characters for readability
                    for i in range(0, len(sequence), 60):
                        f.write(f"{sequence[i:i+60]}\n")
                        
        except (OSError, IOError) as e:
            error_msg = f"Error writing FASTA file {filepath}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileFormatError(error_msg) from e
        except Exception as e:
            error_msg = f"Unexpected error writing FASTA file {filepath}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileFormatError(error_msg) from e
            
        logger.debug(f"Successfully saved FASTA file: {filepath}")
    
    @staticmethod
    def format_excel(df, output_file):
        """
        Save DataFrame to Excel with comprehensive formatting.
        
        Applies professional formatting including headers, color coding,
        column grouping, and proper alignment. Falls back to basic export
        if openpyxl is not available.
        
        Args:
            df: DataFrame with primer results
            output_file: Path to save the formatted Excel file
            
        Returns:
            Path to the saved Excel file
            
        Raises:
            FileFormatError: If Excel file cannot be created
        """
        logger.debug("=== EXCEL FORMATTING DEBUG ===")
        
        # When openpyxl is available, apply full formatting
        try:
            if not HAS_OPENPYXL:
                raise ImportError("openpyxl is not available")
                
            logger.debug(f"Saving formatted Excel file to: {output_file}")
            
            # First, save with pandas to get a basic Excel file
            df.to_excel(output_file, index=False, engine='openpyxl')
            
            # Now open the file for formatting
            workbook = openpyxl.load_workbook(output_file)
            worksheet = workbook.active
            
            # Get the number of rows and columns
            max_row = worksheet.max_row
            max_col = worksheet.max_column
            
            # Create a new row for our custom headers before doing anything else
            worksheet.insert_rows(1)
            max_row += 1
            
            # Create styles
            header_font = Font(bold=True)
            sequence_fill = PatternFill(
                start_color='D9D9D9',  # Light gray
                end_color='D9D9D9',
                fill_type='solid'
            )
            centered_alignment = Alignment(horizontal='center', vertical='center')
            left_alignment = Alignment(horizontal='left', vertical='center')
            
            # Create a dictionary to map column names to their indices
            column_map = {}
            header_texts = []
            
            # First apply basic formatting to all cells before merging
            for col_num in range(1, max_col + 1):
                col_letter = get_column_letter(col_num)
                worksheet.column_dimensions[col_letter].width = 15  # Slightly wider columns
                
                # Set header row formatting in row 2 (original headers)
                cell2 = worksheet.cell(row=2, column=col_num)
                cell2.font = header_font
                cell2.alignment = centered_alignment
                
                # Merge "Gene" header if applicable
                if cell2.value == "Gene":
                    cell1 = worksheet.cell(row=1, column=col_num)
                    cell1.value = "Gene"
                    cell1.font = header_font
                    cell1.alignment = centered_alignment
                    cell1.border = openpyxl.styles.borders.Border()
                    worksheet.merge_cells(start_row=1, start_column=col_num, end_row=2, end_column=col_num)

                # Get the header text and store it
                header_text = cell2.value
                header_texts.append(header_text)
                column_map[header_text] = col_num
                
                # Format data cells
                for row_num in range(3, max_row + 1):
                    cell = worksheet.cell(row=row_num, column=col_num)
                    
                    # Apply sequence fill to Primer F, Primer R, Probe, and Amplicon columns
                    if header_text in ["Primer F", "Primer R", "Probe", "Amplicon"]:
                        cell.fill = sequence_fill
                    
                    # Special handling for "No suitable primers found" cells - left-align these
                    if cell.value == "No suitable primers found":
                        cell.alignment = left_alignment
                    else:
                        # Center alignment for all other cells
                        cell.alignment = centered_alignment
            
            # Freeze panes - freeze first column and first two rows
            worksheet.freeze_panes = 'B3'
            
            # Group the columns for our custom header row
            header_groups = {
                "Gene": [],
                "Forward Primer": ["Primer F", "Tm F", "Penalty F", "Primer F dG", "Primer F BLAST1", "Primer F BLAST2"],
                "Reverse Primer": ["Primer R", "Tm R", "Penalty R", "Primer R dG", "Primer R BLAST1", "Primer R BLAST2", "Pair Penalty"],
                "Probe": ["Probe", "Probe Tm", "Probe Penalty", "Probe dG", "Probe BLAST1", "Probe BLAST2"],
                "Amplicon": ["Amplicon", "Length", "Amplicon GC%", "Amplicon dG"],
                "Location": ["Chromosome", "Location", "Qry Chromosome", "Qry Location"]
            }
            
            # Add the group headers safely by first checking if any columns from each group exist
            for group_name, headers in header_groups.items():
                # Filter to only include headers that actually exist in our DataFrame
                existing_headers = [h for h in headers if h in header_texts]
                
                if not existing_headers:
                    logger.debug(f"Skipping group '{group_name}' as none of its columns exist")
                    continue
                    
                # Find column indices for existing headers
                col_indices = [column_map[h] for h in existing_headers]
                
                if col_indices:
                    start_col = min(col_indices)
                    end_col = max(col_indices)
                    
                    # Add the group header
                    group_cell = worksheet.cell(row=1, column=start_col)
                    group_cell.value = group_name
                    group_cell.font = header_font
                    group_cell.alignment = centered_alignment
                    
                    # Merge if there's more than one column
                    if start_col != end_col:
                        merge_range = f"{get_column_letter(start_col)}1:{get_column_letter(end_col)}1"
                        try:
                            worksheet.merge_cells(merge_range)
                            logger.debug(f"Merged cells for group '{group_name}' ({merge_range})")
                        except Exception as e:
                            # If merge fails, just leave as individual cells
                            logger.warning(f"Could not merge range {merge_range}: {str(e)}")
                            logger.debug(f"Error details: {str(e)}", exc_info=True)
            
            # Save the workbook
            workbook.save(output_file)
            logger.debug(f"Successfully saved formatted Excel file: {output_file}")
            logger.debug("=== END EXCEL FORMATTING DEBUG ===")
            
            return output_file
                
        except ImportError:
            # Expected case when openpyxl is not available
            logger.warning("openpyxl not available, falling back to standard Excel export")
        except Exception as e:
            # If formatting fails, fall back to standard save
            error_msg = f"Error applying Excel formatting: {str(e)}"
            logger.error(error_msg)
            logger.warning("Falling back to standard Excel export")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            
        # Basic fallback save
        try:
            df.to_excel(output_file, index=False)
            logger.info(f"Excel file saved to: {output_file} (without formatting)")
            logger.debug("=== END EXCEL FORMATTING DEBUG ===")
            return output_file
        except Exception as e:
            error_msg = f"Failed to save Excel file {output_file}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END EXCEL FORMATTING DEBUG ===")
            raise FileFormatError(error_msg) from e

    @staticmethod
    def load_sequences_from_table(file_path):
        """
        Load sequences from a CSV or Excel file with improved column detection.
        
        Handles various file formats and column arrangements using intelligent
        analysis to identify sequence name and sequence data columns. Supports
        auto-generation of sequence IDs when names are not provided.
        
        Args:
            file_path: Path to CSV or Excel file
            
        Returns:
            Dictionary of sequence ID to sequence
            
        Raises:
            FileError: If the file doesn't exist
            FileFormatError: If there's an error parsing the file
        """
        logger.debug("=== TABLE LOADING DEBUG ===")
        
        # Validate file exists
        if not os.path.exists(file_path):
            error_msg = f"Table file not found: {file_path}"
            logger.error(error_msg)
            raise FileError(error_msg)
            
        logger.debug(f"Loading sequences from table file: {file_path}")
        
        try:
            # Use SequenceAnalyzer to analyze the file structure
            from ..helpers import SequenceAnalyzer
            analysis = SequenceAnalyzer.analyze_file(file_path)
            
            if "error" in analysis:
                error_msg = f"Error analyzing file: {analysis['error']}"
                logger.error(error_msg)
                raise FileFormatError(error_msg)
                
            # Get recommended columns from analysis
            name_col, seq_col = SequenceAnalyzer.get_recommended_columns(analysis)
            logger.debug(f"Recommended columns: name='{name_col}', sequence='{seq_col}'")
            
            # Determine file type and load with pandas
            try:
                if file_path.endswith('.csv'):
                    df = pd.read_csv(file_path)
                elif file_path.endswith(('.xlsx', '.xls')):
                    df = pd.read_excel(file_path)
                else:
                    error_msg = f"Unsupported file format: {file_path}. Must be CSV or Excel."
                    logger.error(error_msg)
                    raise FileFormatError(error_msg)
            except Exception as e:
                error_msg = f"Error reading file {file_path}: {str(e)}"
                logger.error(error_msg)
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                raise FileFormatError(error_msg) from e
            
            # Remove any empty columns and rows
            df = df.dropna(axis=1, how='all')
            df = df.dropna(axis=0, how='all')
            logger.debug(f"Data shape after cleanup: {df.shape}")
            
            # Convert to dictionary
            sequences = {}
            
            # If no sequence column was found, raise an error
            if not seq_col:
                error_msg = "Could not identify a sequence column in the file"
                logger.error(error_msg)
                raise FileFormatError(error_msg)
                
            # If no name column was found, generate sequential IDs
            if not name_col:
                logger.debug("No sequence name column found, generating sequential IDs")
                for idx, row in df.iterrows():
                    sequence = str(row[seq_col]).strip().upper()
                    
                    # Skip empty sequences
                    if pd.isna(sequence) or not sequence:
                        logger.debug(f"Skipping empty sequence at row {idx}")
                        continue
                    
                    # Generate a sequential ID
                    seq_id = f"Seq_{idx+1}"
                    sequences[seq_id] = sequence
                    logger.debug(f"Added sequence: {seq_id} ({len(sequence)} bp)")
            else:
                # Use the identified name and sequence columns
                logger.debug(f"Using name column '{name_col}' and sequence column '{seq_col}'")
                for idx, row in df.iterrows():
                    seq_id = str(row[name_col]).strip()
                    sequence = str(row[seq_col]).strip().upper()
                    
                    # Skip rows with empty names or sequences
                    if pd.isna(seq_id) or pd.isna(sequence) or not seq_id or not sequence:
                        logger.debug(f"Skipping empty row at index {idx}")
                        continue
                    
                    # Handle duplicate IDs by appending a suffix
                    original_seq_id = seq_id
                    if seq_id in sequences:
                        suffix = 1
                        while f"{seq_id}_{suffix}" in sequences:
                            suffix += 1
                        seq_id = f"{seq_id}_{suffix}"
                        logger.debug(f"Renamed duplicate sequence ID '{original_seq_id}' to '{seq_id}'")
                    
                    sequences[seq_id] = sequence
                    logger.debug(f"Added sequence: {seq_id} ({len(sequence)} bp)")
            
            logger.debug(f"Successfully loaded {len(sequences)} sequences from table file")
            logger.debug("=== END TABLE LOADING DEBUG ===")
            return sequences
            
        except (FileFormatError, FileError):
            # Re-raise specific exceptions without modification
            raise
        except Exception as e:
            error_msg = f"Error parsing sequence table {file_path}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            logger.debug("=== END TABLE LOADING DEBUG ===")
            raise FileFormatError(error_msg) from e

    @classmethod
    def select_sequences_file(cls):
        """
        Prompt the user to select a CSV or Excel file with sequences.
        
        Provides user-friendly file selection specifically for sequence
        table files with appropriate file type filters.
        
        Returns:
            Path to the selected file
            
        Raises:
            FileSelectionError: If file selection fails
        """
        logger.info("\n>>> Please select CSV or Excel file with sequences <<<")
        
        try:
            file_path = cls.select_file(
                "Select CSV or Excel file with sequences", 
                [("CSV Files", "*.csv"), ("Excel Files", "*.xlsx"), ("Excel Files", "*.xls"), ("All Files", "*")]
            )
            logger.debug(f"Selected sequences file: {file_path}")
            return file_path
        except FileSelectionError:
            # Re-raise FileSelectionError without modification
            raise
        except Exception as e:
            error_msg = f"Failed to select sequence file: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileSelectionError(error_msg) from e

    @classmethod
    def select_fasta_file(cls, prompt="Select FASTA file"):
        """
        Prompt the user to select a FASTA file.
        
        Provides user-friendly file selection specifically for FASTA
        files with appropriate file type filters and validation.
        
        Args:
            prompt: Text to display in the file dialog
            
        Returns:
            Path to the selected FASTA file
            
        Raises:
            FileSelectionError: If file selection fails
        """
        try:
            file_path = cls.select_file(
                prompt, 
                [("FASTA Files", "*.fasta"), ("FASTA Files", "*.fa"), ("FASTA Files", "*.fna"), ("All Files", "*")]
            )
            logger.debug(f"Selected FASTA file: {file_path}")
            return file_path
        except FileSelectionError:
            # Re-raise FileSelectionError without modification
            raise
        except Exception as e:
            error_msg = f"Failed to select FASTA file: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileSelectionError(error_msg) from e

    @staticmethod
    def save_results(df, output_dir, input_file, mode='standard', second_fasta=None):
        """
        Save results to an Excel file with correct naming based on input files and mode.
        
        Generates appropriate output filenames based on pipeline mode and input
        files, then saves results with comprehensive formatting.
        
        Args:
            df: DataFrame with primer results
            output_dir: Output directory
            input_file: Path to the input file (FASTA, CSV, MAF, etc.)
            mode: Pipeline mode ('standard', 'direct', or 'alignment')
            second_fasta: Path to second FASTA file for alignment mode
            
        Returns:
            Path to the output file
            
        Raises:
            FileFormatError: If there's an error saving the results
        """
        logger.debug("=== RESULTS SAVING DEBUG ===")
        logger.debug("\nSaving results...")
        
        # Make sure output directory exists
        try:
            os.makedirs(output_dir, exist_ok=True)
        except OSError as e:
            error_msg = f"Failed to create output directory {output_dir}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileFormatError(error_msg) from e
        
        # Set output filename based on mode and input files
        if mode == 'alignment':
            # Alignment mode has specific naming conventions
            if input_file and input_file.endswith('.maf'):
                # Using pre-computed MAF file
                basename = os.path.basename(input_file)
                root, _ = os.path.splitext(basename)
                output_filename = f"Primers_{root}.xlsx"
            elif input_file and second_fasta:
                # Using FASTA alignment with two genomes - use both filenames
                ref_basename = os.path.basename(input_file)
                ref_root, _ = os.path.splitext(ref_basename)
                
                qry_basename = os.path.basename(second_fasta)
                qry_root, _ = os.path.splitext(qry_basename)
                
                output_filename = f"Primers_{ref_root}_vs_{qry_root}.xlsx"
            else:
                # Fallback for alignment mode when files aren't specified properly
                output_filename = f"Primers_Alignment.xlsx"
        elif mode == 'direct':
            # For direct mode, use the basename of the input file
            basename = os.path.basename(input_file)
            root, _ = os.path.splitext(basename)
            output_filename = f"Primers_{root}.xlsx"
        else:
            # Standard mode - use basename of input file
            basename = os.path.basename(input_file)
            root, _ = os.path.splitext(basename)
            output_filename = f"Primers_{root}.xlsx"
        
        # Combine output directory and filename
        output_file = os.path.join(output_dir, output_filename)
        logger.debug(f"Output file will be: {output_file}")
        
        # Define columns based on mode
        if mode == 'direct':
            # Direct mode excludes location columns
            columns = [
                "Gene", "Primer F", "Tm F", "Penalty F", "Primer F dG", "Primer F BLAST1", "Primer F BLAST2",
                "Primer R", "Tm R", "Penalty R", "Primer R dG", "Primer R BLAST1", "Primer R BLAST2",
                "Pair Penalty", "Amplicon", "Length", "Amplicon GC%", "Amplicon dG"
            ]
        else:
            # Standard and alignment modes include location columns
            columns = [
                "Gene", "Primer F", "Tm F", "Penalty F", "Primer F dG", "Primer F BLAST1", "Primer F BLAST2",
                "Primer R", "Tm R", "Penalty R", "Primer R dG", "Primer R BLAST1", "Primer R BLAST2",
                "Pair Penalty", "Amplicon", "Length", "Amplicon GC%", "Amplicon dG",
                "Chromosome", "Location"
            ]
            
            # Add specific columns for alignment mode
            if mode == 'alignment':
                alignment_cols = ["Qry Chromosome", "Qry Location"]
                for col in alignment_cols:
                    if col in df.columns and not df[col].isna().all():
                        columns.append(col)
        
        # Add probe columns if present
        if "Probe" in df.columns:
            probe_cols = [
                "Probe", "Probe Tm", "Probe Penalty", "Probe dG", 
                "Probe BLAST1", "Probe BLAST2"
            ]
            # Insert probe columns after primer columns
            idx = columns.index("Pair Penalty") + 1
            for col in reversed(probe_cols):
                if col in df.columns:
                    columns.insert(idx, col)
            
            logger.debug(f"Added probe columns to output: {', '.join([c for c in probe_cols if c in df.columns])}")
        
        # Ensure all columns in the list exist in the DataFrame
        original_column_count = len(columns)
        columns = [col for col in columns if col in df.columns]
        if len(columns) < original_column_count:
            missing_cols = original_column_count - len(columns)
            logger.debug(f"Filtered out {missing_cols} missing columns from output")
        
        # Reorder the columns
        df = df[columns]
        
        # Save with formatting
        try:
            output_path = FileIO.format_excel(df, output_file)
            logger.debug("=== END RESULTS SAVING DEBUG ===")
            return output_path
        except Exception as e:
            # Fallback if formatting fails
            error_msg = f"Error saving Excel file: {str(e)}"
            logger.error(error_msg)
            logger.warning("Falling back to basic Excel export")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            
            try:
                df.to_excel(output_file, index=False)
                logger.info(f"\nResults saved to: {output_file} (without formatting)")
                logger.debug("=== END RESULTS SAVING DEBUG ===")
                return output_file
            except Exception as ex:
                error_msg = f"Failed to save results to {output_file}: {str(ex)}"
                logger.error(error_msg)
                logger.debug(f"Error details: {str(ex)}", exc_info=True)
                logger.debug("=== END RESULTS SAVING DEBUG ===")
                raise FileFormatError(error_msg) from ex

    @staticmethod
    def setup_output_directories(args, reference_file=None, mode='standard'):
        """
        Set up output directory and temporary directory for the pipeline.
        
        Creates appropriate directory structure based on pipeline mode and
        input files, with intelligent path determination for different scenarios.
        
        Args:
            args: Command line arguments
            reference_file: Path to reference file
            mode: Pipeline mode ('standard', 'direct', or 'alignment')
            
        Returns:
            Tuple of (output_dir, temp_dir) - Paths to output and temporary directories
            
        Raises:
            FileError: If directory creation fails
        """
        logger.debug("=== OUTPUT DIRECTORY SETUP DEBUG ===")
        
        # Determine output directory
        if args.output:
            output_dir = args.output
            logger.debug(f"Using user-specified output directory: {output_dir}")
        elif reference_file:
            if mode == 'alignment' and hasattr(args, 'maf_file') and args.maf_file:
                # For MAF files, use the parent directory of the "Alignments" folder
                maf_dir = os.path.dirname(os.path.abspath(args.maf_file))
                if os.path.basename(maf_dir) == "Alignments":
                    # If it's in an Alignments folder, go up one level
                    output_dir = os.path.join(os.path.dirname(maf_dir), "Primers")
                else:
                    # Otherwise use the parent directory of the MAF file
                    output_dir = os.path.join(os.path.dirname(maf_dir), "Primers")
                logger.debug(f"MAF-based output directory: {output_dir}")
            else:
                # Use the directory of the reference file
                input_dir = os.path.dirname(os.path.abspath(reference_file))
                output_dir = os.path.join(input_dir, "Primers")
                logger.debug(f"Reference file-based output directory: {output_dir}")
        else:
            # Fallback to current directory if no reference file
            output_dir = os.path.join(os.getcwd(), "Primers")
            logger.debug(f"Fallback output directory: {output_dir}")
        
        # Create output directory
        try:
            os.makedirs(output_dir, exist_ok=True)
            logger.debug(f"Created output directory: {output_dir}")
        except OSError as e:
            error_msg = f"Failed to create output directory {output_dir}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(error_msg) from e
        
        # Create temporary directory
        try:
            temp_dir = tempfile.mkdtemp(prefix="ddprimer_temp_", dir=output_dir)
            logger.debug(f"Created temporary directory: {temp_dir}")
        except OSError as e:
            error_msg = f"Failed to create temporary directory in {output_dir}: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(error_msg) from e
        
        logger.debug("=== END OUTPUT DIRECTORY SETUP DEBUG ===")
        return output_dir, temp_dir

    @staticmethod
    def cleanup_temp_directory(temp_dir):
        """
        Clean up temporary directory safely.
        
        Removes temporary directory and all contents with comprehensive
        error handling and logging.
        
        Args:
            temp_dir: Path to temporary directory
            
        Returns:
            True if cleanup was successful, False otherwise
        """
        if not temp_dir:
            logger.debug("No temporary directory to clean up")
            return True
            
        try:
            if os.path.exists(temp_dir):
                logger.debug(f"Cleaning up temporary directory: {temp_dir}")
                shutil.rmtree(temp_dir)
                logger.debug("Temporary directory cleanup completed successfully")
                return True
            else:
                logger.debug(f"Temporary directory does not exist: {temp_dir}")
                return True
        except OSError as e:
            logger.warning(f"Error cleaning up temporary files in {temp_dir}: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return False
        except Exception as e:
            logger.warning(f"Unexpected error cleaning up temporary files: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return False
    

class TempDirectoryManager:
    """
    Context manager for temporary directory creation and cleanup.
    
    Provides safe temporary directory management with automatic cleanup
    and comprehensive error handling for pipeline operations.
    
    Attributes:
        temp_dir: Path to the created temporary directory
        base_dir: Base directory for temporary directory creation
        
    Example:
        >>> with TempDirectoryManager() as temp_dir:
        ...     # Use temp_dir for operations
        ...     process_files(temp_dir)
        # temp_dir is automatically cleaned up
    """
    
    def __init__(self, base_dir=None):
        """
        Initialize the temporary directory manager.
        
        Args:
            base_dir: Base directory to create the temp directory in
        """
        self.temp_dir = None
        self.base_dir = base_dir
        
    def __enter__(self):
        """
        Create and return the temporary directory path.
        
        Returns:
            Path to the temporary directory
            
        Raises:
            FileError: If temporary directory creation fails
        """
        try:
            self.temp_dir = tempfile.mkdtemp(prefix="ddprimer_temp_", dir=self.base_dir)
            logger.debug(f"Created temporary directory: {self.temp_dir}")
            return self.temp_dir
        except OSError as e:
            error_msg = f"Failed to create temporary directory: {str(e)}"
            logger.error(error_msg)
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileError(error_msg) from e
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Clean up the temporary directory.
        
        Args:
            exc_type: Exception type if an exception was raised
            exc_val: Exception value if an exception was raised
            exc_tb: Exception traceback if an exception was raised
        """
        if self.temp_dir and os.path.exists(self.temp_dir):
            try:
                logger.debug(f"Cleaning up temporary directory: {self.temp_dir}")
                shutil.rmtree(self.temp_dir)
                logger.debug("Temporary directory cleanup completed")
            except OSError as e:
                logger.warning(f"Error cleaning up temporary files in {self.temp_dir}: {str(e)}")
                logger.debug(f"Error details: {str(e)}", exc_info=True)
            except Exception as e:
                logger.warning(f"Unexpected error during cleanup: {str(e)}")
                logger.debug(f"Error details: {str(e)}", exc_info=True)