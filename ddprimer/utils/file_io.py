#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File I/O module for ddPrimer pipeline.

This module provides consistent interfaces for file operations including:
- File selection
- File reading/writing
- Format conversions
"""

import os
import sys
import platform
import logging
import pandas as pd
import contextlib
import tempfile
import shutil

from ..config import Config
from ..config.exceptions import FileSelectionError, FileFormatError


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
    """Temporarily redirect C-level stderr to /dev/null (macOS IMK chatter)."""
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
    """Handles all file I/O operations with consistent interfaces."""
    
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
        """
        if HAS_WX and cls._wx_app is None:
            with _silence_cocoa_stderr():
                cls._wx_app = wx.App(False)
                logger = logging.getLogger("ddPrimer")
                logger.debug("wxPython app initialized")
    
    @classmethod
    def hide_app(cls):
        """
        Hide the wxPython app without destroying it.
        This is safer on macOS and avoids segmentation faults.
        """
        if cls.is_macos and cls.has_pyobjc and cls._wx_app is not None:
            try:
                import AppKit
                # Hide from Dock
                NSApplication = AppKit.NSApplication.sharedApplication()
                NSApplication.setActivationPolicy_(AppKit.NSApplicationActivationPolicyAccessory)
                logger = logging.getLogger("ddPrimer")
                logger.debug("App hidden from macOS dock")
            except Exception as e:
                logger = logging.getLogger("ddPrimer")
                logger.debug(f"Error hiding app: {e}")

    @classmethod
    def mark_selection_complete(cls):
        """
        Mark that all file selections are complete and hide the app from the dock.
        """
        logger = logging.getLogger("ddPrimer")
        logger.debug("File selection process marked as complete")
        cls.hide_app()

    @classmethod
    def get_last_directory(cls):
        """
        Get the last directory used, loading from persistent storage if needed.
        
        Returns:
            str: Path to the last directory used
        """
        if cls._last_directory is None:
            # Try to load the last directory from a config file
            logger = logging.getLogger("ddPrimer")
            config_dir = os.path.join(os.path.expanduser("~"), ".ddPrimer")
            os.makedirs(config_dir, exist_ok=True)
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
                except Exception as e:
                    logger.debug(f"Error reading last directory config: {e}")
                    cls._last_directory = os.path.expanduser("~")
            else:
                cls._last_directory = os.path.expanduser("~")
        
        return cls._last_directory
    
    @classmethod
    def save_last_directory(cls, directory):
        """
        Save the last directory to persistent storage.
        
        Args:
            directory (str): Path to save
        """
        if directory and os.path.isdir(directory):
            cls._last_directory = directory
            logger = logging.getLogger("ddPrimer")
            
            # Save to config file
            config_dir = os.path.join(os.path.expanduser("~"), ".ddPrimer")
            os.makedirs(config_dir, exist_ok=True)
            config_file = os.path.join(config_dir, "last_directory.txt")
            
            try:
                with open(config_file, 'w') as f:
                    f.write(directory)
                logger.debug(f"Saved last directory to config: {directory}")
            except Exception as e:
                logger.debug(f"Error saving last directory config: {e}")

    @classmethod
    def normalize_filetypes(cls, filetypes):
        """
        Normalize file types for different platforms.
        
        Args:
            filetypes (list): List of (description, extension) tuples
            
        Returns:
            list: Normalized file types list
        """
        logger = logging.getLogger("ddPrimer")
        
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
                    
        logger.debug(f"Normalized {len(filetypes) if filetypes else 0} file types for platform {platform.system()}")
        return normalized
        
    @classmethod
    def select_file(cls, prompt, filetypes):
        """
        Show a file dialog or prompt for a file path if CLI mode is enabled.
        
        Args:
            prompt (str): Text to display in the file dialog
            filetypes (list): List of file type tuples for the dialog
            
        Returns:
            str: Selected file path
            
        Raises:
            FileSelectionError: If file selection fails
        """
        logger = logging.getLogger("ddPrimer")
        
        # Get the last directory from our persistent storage
        last_directory = cls.get_last_directory()
        logger.debug(f"Starting file selection from directory: {last_directory}")

        if cls.use_cli:
            logger.debug(f"{prompt} (CLI mode)")
            file_path = input(f"{prompt}: ").strip()
            if not file_path:
                logger.error("No file path provided in CLI mode")
                raise FileSelectionError("No file path provided")
                
            # Validate the path exists
            if not os.path.exists(file_path):
                logger.error(f"File not found: {file_path}")
                raise FileSelectionError(f"File not found: {file_path}")
                
            return file_path

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
                    logger.error("No file was selected in the dialog")
                    raise FileSelectionError("No file was selected")

                # Update the last directory
                cls.save_last_directory(os.path.dirname(file_path))
                logger.debug(f"Selected file: {file_path}")
                return file_path

            except Exception as e:
                logger.error(f"wxPython file selection failed: {e}")
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
        
        Args:
            prompt (str): Text to display in the file dialog
            filetypes (list): List of file type tuples for the dialog
            
        Returns:
            list: List of selected file paths
            
        Raises:
            FileSelectionError: If file selection fails
        """
        logger = logging.getLogger("ddPrimer")
        
        # Get the last directory from our persistent storage
        last_directory = cls.get_last_directory()
        logger.debug(f"Starting multiple file selection from directory: {last_directory}")

        if cls.use_cli:
            logger.info(f"{prompt} (CLI mode)")
            paths = input(f"{prompt} (paths separated by spaces): ").strip()
            file_paths = paths.split() if paths else []
            if not file_paths:
                logger.error("No file paths provided in CLI mode")
                raise FileSelectionError("No file paths provided")
                
            # Validate the paths exist
            for path in file_paths:
                if not os.path.exists(path):
                    logger.error(f"File not found: {path}")
                    raise FileSelectionError(f"File not found: {path}")
                    
            return file_paths

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
                    logger.error("No files were selected in the dialog")
                    raise FileSelectionError("No files were selected")

                # Update the last directory
                cls.save_last_directory(os.path.dirname(file_paths[0]))
                logger.debug(f"Selected {len(file_paths)} file(s)")
                return file_paths

            except Exception as e:
                logger.error(f"wxPython multiple file selection failed: {e}")
                logger.debug(f"Error details: {str(e)}", exc_info=True)
                # Fall back to CLI mode if wxPython fails
                cls.use_cli = True

        # If wxPython isn't available or also failed, use CLI
        logger.warning("Falling back to CLI mode for multiple file selection")
        return cls.select_files(prompt, filetypes)
        
    @staticmethod
    def load_fasta(filepath):
        """
        Load sequences from a FASTA file into a dict: {header_without_gt: sequence}.
        Optimized for memory efficiency.
        
        Args:
            filepath (str): Path to the FASTA file
            
        Returns:
            dict: Dictionary mapping sequence headers to sequences
            
        Raises:
            FileNotFoundError: If the FASTA file doesn't exist
            FileFormatError: If there's an error parsing the FASTA file
        """
        logger = logging.getLogger("ddPrimer")
        sequences = {}
        name = None
        seq_chunks = []
        
        logger.debug(f"Loading FASTA file: {filepath}")
        
        # Validate file exists
        if not os.path.exists(filepath):
            logger.error(f"FASTA file not found: {filepath}")
            raise FileNotFoundError(f"FASTA file not found: {filepath}")
            
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if name:
                            sequences[name] = "".join(seq_chunks).upper()
                        name = line[1:].split()[0]
                        seq_chunks = []
                    else:
                        seq_chunks.append(line)
                if name:
                    sequences[name] = "".join(seq_chunks).upper()
        except Exception as e:
            logger.error(f"Error reading FASTA file {os.path.abspath(filepath)}: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileFormatError(f"Error parsing FASTA file: {str(e)}")
        
        logger.debug(f"Loaded {len(sequences)} sequences from FASTA file")
        return sequences
    
    @staticmethod
    def save_fasta(sequences, filepath):
        """
        Save sequences to a FASTA file.
        
        Args:
            sequences (dict): Dictionary of sequence ID to sequence
            filepath (str): Path to save the FASTA file
            
        Raises:
            FileFormatError: If there's an error writing the FASTA file
        """
        logger = logging.getLogger("ddPrimer")
        logger.info(f"Saving {len(sequences)} sequences to FASTA file: {filepath}")
        
        try:
            with open(filepath, 'w') as f:
                for seq_id, sequence in sequences.items():
                    f.write(f">{seq_id}\n")
                    # Write sequence in chunks of 60 characters for readability
                    for i in range(0, len(sequence), 60):
                        f.write(f"{sequence[i:i+60]}\n")
        except Exception as e:
            logger.error(f"Error writing FASTA file {filepath}: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileFormatError(f"Error writing FASTA file: {str(e)}")
            
        logger.debug(f"Successfully saved FASTA file: {filepath}")
    
    @staticmethod
    def format_excel(df, output_file):
        """
        Save DataFrame to Excel with specific formatting.
        If openpyxl is not available, falls back to standard pandas Excel export.
        
        Args:
            df (pandas.DataFrame): DataFrame with primer results
            output_file (str): Path to save the formatted Excel file
            
        Returns:
            str: Path to the saved Excel file
        """
        logger = logging.getLogger("ddPrimer")
        
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
                            logger.debug(f"Could not merge range {merge_range}: {str(e)}")
            
            # Save the workbook
            workbook.save(output_file)
            logger.debug(f"Successfully saved formatted Excel file: {output_file}")
            
            return output_file
                
        except Exception as e:
            # If formatting fails, fall back to standard save
            logger.error(f"Error applying Excel formatting: {str(e)}")
            logger.warning("Falling back to standard Excel export")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            
            # Basic fallback save
            df.to_excel(output_file, index=False)
            logger.info(f"Excel file saved to: {output_file} (without formatting)")
            
            return output_file

        
    @staticmethod
    def load_sequences_from_table(file_path):
        """
        Load sequences from a CSV or Excel file with improved column detection.
        
        This method can handle various file formats:
        - Files with explicit sequence name and sequence columns
        - Files with only sequence data (auto-generated IDs)
        - Files with multiple columns where sequence columns aren't adjacent to name columns
        
        Args:
            file_path (str): Path to CSV or Excel file
            
        Returns:
            dict: Dictionary of sequence ID to sequence
            
        Raises:
            FileNotFoundError: If the file doesn't exist
            FileFormatError: If there's an error parsing the file
        """
        logger = logging.getLogger("ddPrimer")
        
        # Validate file exists
        if not os.path.exists(file_path):
            logger.error(f"Table file not found: {file_path}")
            raise FileNotFoundError(f"Table file not found: {file_path}")
            
        logger.debug(f"Loading sequences from table file: {file_path}")
        
        try:
            # Use SequenceAnalyzer to analyze the file structure
            from ..helpers.sequence_analyzer import SequenceAnalyzer
            analysis = SequenceAnalyzer.analyze_file(file_path)
            
            if "error" in analysis:
                raise FileFormatError(f"Error analyzing file: {analysis['error']}")
                
            # Get recommended columns from analysis
            name_col, seq_col = SequenceAnalyzer.get_recommended_columns(analysis)
            
            # Determine file type and load with pandas
            if file_path.endswith('.csv'):
                df = pd.read_csv(file_path)
            elif file_path.endswith(('.xlsx', '.xls')):
                df = pd.read_excel(file_path)
            else:
                logger.error(f"Unsupported file format: {file_path}")
                raise FileFormatError(f"Unsupported file format: {file_path}. Must be CSV or Excel.")
            
            # Remove any empty columns and rows
            df = df.dropna(axis=1, how='all')
            df = df.dropna(axis=0, how='all')
            
            # Convert to dictionary
            sequences = {}
            
            # If no sequence column was found, raise an error
            if not seq_col:
                logger.error("Could not identify a sequence column in the file")
                raise FileFormatError("Could not identify a sequence column in the file")
                
            # If no name column was found, generate sequential IDs
            if not name_col:
                logger.debug("No sequence name column found, generating sequential IDs")
                for idx, row in df.iterrows():
                    sequence = str(row[seq_col]).strip().upper()
                    
                    # Skip empty sequences
                    if pd.isna(sequence) or not sequence:
                        continue
                    
                    # Generate a sequential ID
                    seq_id = f"Seq_{idx+1}"
                    sequences[seq_id] = sequence
            else:
                # Use the identified name and sequence columns
                for idx, row in df.iterrows():
                    seq_id = str(row[name_col]).strip()
                    sequence = str(row[seq_col]).strip().upper()
                    
                    # Skip rows with empty names or sequences
                    if pd.isna(seq_id) or pd.isna(sequence) or not seq_id or not sequence:
                        continue
                    
                    # Handle duplicate IDs by appending a suffix
                    if seq_id in sequences:
                        suffix = 1
                        while f"{seq_id}_{suffix}" in sequences:
                            suffix += 1
                        seq_id = f"{seq_id}_{suffix}"
                        logger.debug(f"Renamed duplicate sequence ID to '{seq_id}'")
                    
                    sequences[seq_id] = sequence
            
            logger.debug(f"Loaded {len(sequences)} sequences from table file")
            return sequences
            
        except Exception as e:
            if isinstance(e, FileFormatError) or isinstance(e, FileNotFoundError):
                raise e
            logger.error(f"Error loading sequences from table file: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileFormatError(f"Error parsing sequence table: {str(e)}")

    @classmethod
    def select_sequences_file(cls):
        """
        Prompt the user to select a CSV or Excel file with sequences.
        
        Returns:
            str: Path to the selected file
            
        Raises:
            FileSelectionError: If file selection fails
        """
        logger = logging.getLogger("ddPrimer")
        
        logger.info("\n>>> Please select CSV or Excel file with sequences <<<")
        
        try:
            file_path = cls.select_file(
                "Select CSV or Excel file with sequences", 
                [("CSV Files", "*.csv"), ("Excel Files", "*.xlsx"), ("Excel Files", "*.xls"), ("All Files", "*")]
            )
            logger.debug(f"Selected sequences file: {file_path}")
            return file_path
        except Exception as e:
            logger.error(f"Error selecting sequence file: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileSelectionError(f"Failed to select sequence file: {str(e)}")

    @classmethod
    def select_fasta_file(cls, prompt="Select FASTA file"):
        """
        Prompt the user to select a FASTA file.
        
        Args:
            prompt (str): Text to display in the file dialog
            
        Returns:
            str: Path to the selected FASTA file
            
        Raises:
            FileSelectionError: If file selection fails
        """
        logger = logging.getLogger("ddPrimer")

        
        try:
            file_path = cls.select_file(
                prompt, 
                [("FASTA Files", "*.fasta"), ("FASTA Files", "*.fa"), ("FASTA Files", "*.fna"), ("All Files", "*")]
            )
            logger.debug(f"Selected FASTA file: {file_path}")
            return file_path
        except Exception as e:
            logger.error(f"Error selecting FASTA file: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            raise FileSelectionError(f"Failed to select FASTA file: {str(e)}")

    @staticmethod
    def save_results(df, output_dir, input_file, mode='standard', second_fasta=None):
        """
        Save results to an Excel file with correct naming based on input files and mode.
        
        Args:
            df (pandas.DataFrame): DataFrame with primer results
            output_dir (str): Output directory
            input_file (str): Path to the input file (FASTA, CSV, MAF, etc.)
            mode (str): Pipeline mode ('standard', 'direct', or 'alignment')
            second_fasta (str, optional): Path to second FASTA file for alignment mode
            
        Returns:
            str: Path to the output file
            
        Raises:
            FileFormatError: If there's an error saving the results
        """
        logger = logging.getLogger("ddPrimer")
            
        logger.debug("\nSaving results...")
        
        # Make sure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
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
            
            # Add cross-species specific columns for alignment mode
            if mode == 'alignment':
                cross_species_cols = ["Qry Chromosome", "Qry Location"]
                for col in cross_species_cols:
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
        columns = [col for col in columns if col in df.columns]
        
        # Reorder the columns
        df = df[columns]
        
        # Save with formatting
        try:
            output_path = FileIO.format_excel(df, output_file)
            return output_path
        except Exception as e:
            # Fallback if formatting fails
            logger.error(f"Error saving Excel file: {str(e)}")
            logger.warning("Falling back to basic Excel export")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            
            try:
                df.to_excel(output_file, index=False)
                logger.info(f"\nResults saved to: {output_file} (without formatting)")
                return output_file
            except Exception as ex:
                logger.error(f"Failed to save results: {str(ex)}")
                logger.debug(f"Error details: {str(ex)}", exc_info=True)
                raise FileFormatError(f"Failed to save results: {str(ex)}")

    @staticmethod
    def setup_output_directories(args, reference_file=None, mode='standard'):
        """
        Set up output directory and temporary directory for the pipeline.
        
        Args:
            args (argparse.Namespace): Command line arguments
            reference_file (str, optional): Path to reference file
            mode (str): Pipeline mode ('standard', 'direct', or 'alignment')
            
        Returns:
            tuple: (output_dir, temp_dir) - Paths to output and temporary directories
        """
        logger = logging.getLogger("ddPrimer")
        
        # Determine output directory
        if args.output:
            output_dir = args.output
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
            else:
                # Use the directory of the reference file
                input_dir = os.path.dirname(os.path.abspath(reference_file))
                output_dir = os.path.join(input_dir, "Primers")
        else:
            # Fallback to current directory if no reference file
            output_dir = os.path.join(os.getcwd(), "Primers")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        logger.debug(f"Created output directory: {output_dir}")
        
        # Create temporary directory
        temp_dir = tempfile.mkdtemp(prefix="ddprimer_temp_", dir=output_dir)
        logger.debug(f"Created temporary directory: {temp_dir}")
        
        return output_dir, temp_dir

    @staticmethod
    def cleanup_temp_directory(temp_dir):
        """
        Clean up temporary directory safely.
        
        Args:
            temp_dir (str): Path to temporary directory
            
        Returns:
            bool: True if cleanup was successful, False otherwise
        """
        logger = logging.getLogger("ddPrimer")
        
        if not temp_dir:
            return True
            
        try:
            if os.path.exists(temp_dir):
                logger.debug(f"Cleaning up temporary directory: {temp_dir}")
                shutil.rmtree(temp_dir)
                return True
        except Exception as e:
            logger.warning(f"Error cleaning up temporary files: {str(e)}")
            logger.debug(f"Error details: {str(e)}", exc_info=True)
            return False
        
        return True
    
class TempDirectoryManager:
    """Context manager for temporary directory creation and cleanup."""
    
    def __init__(self, base_dir=None):
        """
        Initialize the temporary directory manager.
        
        Args:
            base_dir (str, optional): Base directory to create the temp directory in
        """
        self.temp_dir = None
        self.base_dir = base_dir
        self.logger = logging.getLogger("ddPrimer")
        
    def __enter__(self):
        """
        Create and return the temporary directory path.
        
        Returns:
            str: Path to the temporary directory
        """
        self.temp_dir = tempfile.mkdtemp(prefix="ddprimer_temp_", dir=self.base_dir)
        self.logger.debug(f"Created temporary directory: {self.temp_dir}")
        return self.temp_dir
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Clean up the temporary directory.
        
        Args:
            exc_type: Exception type if an exception was raised
            exc_val: Exception value if an exception was raised
            exc_tb: Exception traceback if an exception was raised
        """
        try:
            if self.temp_dir and os.path.exists(self.temp_dir):
                self.logger.debug(f"Cleaning up temporary directory: {self.temp_dir}")
                shutil.rmtree(self.temp_dir)
        except Exception as e:
            self.logger.warning(f"Error cleaning up temporary files: {str(e)}")
            self.logger.debug(f"Error details: {str(e)}", exc_info=True)