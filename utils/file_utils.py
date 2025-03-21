#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File handling utilities for ddPrimer pipeline.
"""

import os
import sys
import platform
import tkinter as tk
from tkinter import filedialog
from ..config import Config


class FileUtils:
    """Handles file operations including UI dialogs and file parsing."""
    
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

    @classmethod
    def get_last_directory(cls):
        """
        Get the last directory used, loading from persistent storage if needed.
        
        Returns:
            str: Path to the last directory used
        """
        if cls._last_directory is None:
            # Try to load the last directory from a config file
            config_dir = os.path.join(os.path.expanduser("~"), ".ddPrimer")
            os.makedirs(config_dir, exist_ok=True)
            config_file = os.path.join(config_dir, "last_directory.txt")
            
            if os.path.exists(config_file):
                try:
                    with open(config_file, 'r') as f:
                        last_dir = f.read().strip()
                    if os.path.isdir(last_dir):
                        cls._last_directory = last_dir
                        Config.debug(f"Loaded last directory from config: {last_dir}")
                    else:
                        cls._last_directory = os.path.expanduser("~")
                except Exception as e:
                    Config.debug(f"Error reading last directory config: {e}")
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
            
            # Save to config file
            config_dir = os.path.join(os.path.expanduser("~"), ".ddPrimer")
            os.makedirs(config_dir, exist_ok=True)
            config_file = os.path.join(config_dir, "last_directory.txt")
            
            try:
                with open(config_file, 'w') as f:
                    f.write(directory)
                Config.debug(f"Saved last directory to config: {directory}")
            except Exception as e:
                Config.debug(f"Error saving last directory config: {e}")

    @staticmethod
    def _show_macos_file_dialog(prompt, filetypes, multiple=False):
        """
        Show a file dialog using tkinter (fallback for macOS).
        
        Args:
            prompt (str): Title for the dialog
            filetypes (list): List of file type tuples for the dialog
            multiple (bool): Allow multiple file selection
            
        Returns:
            str or list: Selected file path(s) or None if cancelled
        """
        try:
            import tkinter as tk
            from tkinter import filedialog

            root = tk.Tk()
            root.withdraw()

            # Redirect stderr to suppress warnings
            original_stderr = sys.stderr
            sys.stderr = open(os.devnull, "w")

            if multiple:
                file_paths = filedialog.askopenfilenames(title=prompt, filetypes=filetypes, parent=root)
            else:
                file_paths = filedialog.askopenfilename(title=prompt, filetypes=filetypes, parent=root)

            sys.stderr.close()
            sys.stderr = original_stderr
            root.destroy()

            if not file_paths:
                return None if not multiple else []

            return file_paths
        except Exception as e:
            Config.debug(f"Error in fallback tkinter dialog: {e}")
            return None if not multiple else []

    @classmethod
    def get_file(cls, prompt, filetypes):
        """
        Show a file dialog or prompt for a file path if CLI mode is enabled.
        
        Args:
            prompt (str): Text to display in the file dialog
            filetypes (list): List of file type tuples for the dialog
            
        Returns:
            str: Selected file path
        """
        # Get the last directory from our persistent storage
        last_directory = cls.get_last_directory()
        Config.debug(f"Starting with last directory: {last_directory}")
        
        if cls.use_cli:
            file_path = input(f"{prompt} (enter file path manually): ").strip()
            if not file_path:
                print(f"{prompt} - No file selected. Exiting.")
                sys.exit(1)
            return file_path
        
        # Fall back to tkinter implementation
        try:
            # Ensure we're starting fresh - kill any existing tk instances
            try:
                if hasattr(tk, '_default_root') and tk._default_root:
                    tk._default_root.destroy()
                    tk._default_root = None
            except:
                pass  # Ignore if no existing root
            
            # Create completely new tk root
            root = tk.Tk()
            root.withdraw()  # Hide the main window
            
            # Let the window manager position the dialog
            # This is more reliable than trying to position it ourselves
            
            # On macOS, we must use a simpler filetype pattern
            valid_filetypes = [("All Files", "*")]
            
            # Add more specific filetypes if not on macOS or if they're simple enough
            if not cls.is_macos:
                for desc, ext in filetypes:
                    if ext and isinstance(ext, str) and ext.strip() and ext != "*":
                        valid_filetypes.insert(0, (desc, ext))
            
            # Make sure dialog doesn't output errors to terminal
            original_stderr = sys.stderr
            sys.stderr = open(os.devnull, "w")
            
            # Log the initial directory we're using
            Config.debug(f"Opening dialog with initial directory: {last_directory}")
            
            try:
                # Show the dialog
                file_path = filedialog.askopenfilename(
                    title=prompt,
                    filetypes=valid_filetypes,
                    initialdir=last_directory
                )
            finally:
                # Restore stderr
                sys.stderr.close()
                sys.stderr = original_stderr
            
            # Force proper cleanup
            root.destroy()
            
            if not file_path:
                print(f"No file was selected in the dialog. Exiting.")
                sys.exit(1)
            
            # Update the last directory for next time
            if os.path.isfile(file_path):
                new_directory = os.path.dirname(file_path)
                cls.save_last_directory(new_directory)
                Config.debug(f"Updated last directory to: {new_directory}")
            
            print(f"File selected: {file_path}")
            return file_path
            
        except Exception as e:
            print(f"GUI initialization failed ({e}), falling back to CLI mode.")
            cls.use_cli = True
            return cls.get_file(prompt, filetypes)

    @classmethod
    def get_files(cls, prompt, filetypes):
        """
        Show a file dialog for multiple files or prompt if CLI mode is enabled.
        
        Args:
            prompt (str): Text to display in the file dialog
            filetypes (list): List of file type tuples for the dialog
            
        Returns:
            list: List of selected file paths
        """
        # Get the last directory from our persistent storage
        last_directory = cls.get_last_directory()
        Config.debug(f"Starting with last directory: {last_directory}")
            
        if cls.use_cli:
            paths = input(f"{prompt} (enter file paths separated by spaces): ").strip()
            file_paths = paths.split() if paths else []
            if not file_paths:
                print(f"{prompt} - No files selected. Exiting.")
                sys.exit(1)
            return file_paths
        
        # Fall back to tkinter implementation
        try:
            # Ensure we're starting fresh - kill any existing tk instances
            try:
                if hasattr(tk, '_default_root') and tk._default_root:
                    tk._default_root.destroy()
                    tk._default_root = None
            except:
                pass  # Ignore if no existing root
            
            # Create completely new tk root
            root = tk.Tk()
            root.withdraw()  # Hide the main window
            
            # Let the window manager position the dialog
            # This is more reliable than trying to position it ourselves
            
            # On macOS, we must use a simpler filetype pattern
            valid_filetypes = [("All Files", "*")]
            
            # Add more specific filetypes if not on macOS or if they're simple enough
            if not cls.is_macos:
                for desc, ext in filetypes:
                    if ext and isinstance(ext, str) and ext.strip() and ext != "*":
                        valid_filetypes.insert(0, (desc, ext))
            
            # Make sure dialog doesn't output errors to terminal
            original_stderr = sys.stderr
            sys.stderr = open(os.devnull, "w")
            
            # Log the initial directory we're using
            Config.debug(f"Opening dialog with initial directory: {last_directory}")
            
            try:
                # Show the dialog
                file_paths = filedialog.askopenfilenames(
                    title=prompt,
                    filetypes=valid_filetypes,
                    initialdir=last_directory
                )
            finally:
                # Restore stderr
                sys.stderr.close()
                sys.stderr = original_stderr
            
            # Force proper cleanup
            root.destroy()
            
            if not file_paths:
                print(f"No files were selected in the dialog. Exiting.")
                sys.exit(1)
            
            # Update the last directory for next time
            if file_paths and os.path.isfile(file_paths[0]):
                new_directory = os.path.dirname(file_paths[0])
                cls.save_last_directory(new_directory)
                Config.debug(f"Updated last directory to: {new_directory}")
            
            print(f"Selected {len(file_paths)} file(s)")
            for path in file_paths:
                print(f" - {path}")
                
            return file_paths
            
        except Exception as e:
            print(f"GUI initialization failed ({e}), falling back to CLI mode.")
            cls.use_cli = True
            return cls.get_files(prompt, filetypes)
    
    @staticmethod
    def load_fasta(filepath):
        """
        Load sequences from a FASTA file into a dict: {header_without_gt: sequence}.
        Optimized for memory efficiency.
        
        Args:
            filepath (str): Path to the FASTA file
            
        Returns:
            dict: Dictionary mapping sequence headers to sequences
        """
        sequences = {}
        name = None
        seq_chunks = []
        
        Config.debug(f"Loading FASTA file: {filepath}")
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
            print(f"Error reading FASTA file {os.path.abspath(filepath)}: {e}")
            sys.exit(1)
        
        Config.debug(f"Loaded {len(sequences)} sequences from FASTA file")
        return sequences