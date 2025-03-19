#!/usr/bin/env python3
"""
ddPrimer Setup Tool
------------------
This script automates the installation and configuration of dependencies 
required for the ddPrimer tool.
"""

import os
import sys
import subprocess
import shutil
import platform
import urllib.request
import logging
from pathlib import Path
import re
import tempfile
import tarfile
import zipfile

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(message)s'
)
logger = logging.getLogger('ddprimer_setup')

# Required Python Packages
PYTHON_PACKAGES = ["pandas", "tqdm", "nupack"]

# System Packages mapping: command_name -> package_name
SYSTEM_PACKAGES = {
    "Linux": {
        "primer3_core": "primer3",
        "blastn": "ncbi-blast+",
        "lastz": "lastz",
        "samtools": "samtools"
    },
    "Darwin": {  # macOS
        "primer3_core": "primer3",
        "blastn": "blast",
        "lastz": "lastz",
        "samtools": "samtools"
    },
    "Windows": {
        # Windows packages are handled differently
    }
}

# PLastZ paths
PLASTZ_DIR = "PLastZ"
PLASTZ_SCRIPT = "PLastZ.py"
SUPPORT_DIR = "/Library/Application Support/ddPrimer"

class DdPrimerSetup:
    """Main class for setting up ddPrimer dependencies."""
    
    def __init__(self):
        """Initialize setup environment variables."""
        self.os_type = platform.system()
        self.blast_db_path = None
        self.multiz_path = None
        self.plastz_path = None
        self.install_dir = os.path.expanduser("~/.ddprimer")
        
        # Properly detect headless mode - set to True if environment variable exists or no display available
        self.headless = os.environ.get('DDPRIMER_HEADLESS') == '1' or not self._is_gui_available()
        
        # Set default paths for headless mode
        if self.headless:
            # Check for environment variables that specify paths
            self.default_multiz = os.environ.get('DDPRIMER_MULTIZ_PATH')
            self.default_plastz = os.environ.get('DDPRIMER_PLASTZ_PATH')
            self.default_fasta = os.environ.get('DDPRIMER_FASTA_PATH')
            
            logger.info("Running in headless mode. GUI dialogs will be disabled.")
    
    def _is_gui_available(self):
        """Check if a GUI is available."""
        # Check for common headless indicators
        if 'DISPLAY' not in os.environ and platform.system() != 'Windows':
            return False
        
        # On some systems, DISPLAY might be set but still be headless (like some HPC nodes)
        if os.environ.get('SSH_TTY') and not os.environ.get('DDPRIMER_GUI_FORCE', False):
            return False
            
        # Try to import tkinter without creating a window
        try:
            import tkinter as tk
            test = tk.Tk()
            test.withdraw()
            test.destroy()
            return True
        except (ImportError, tk.TclError):
            return False
    
    def run(self):
        """Run the complete setup process."""
        logger.info("=== Setting up ddPrimer Dependencies ===")
        logger.info(f"Running in {'headless' if self.headless else 'GUI'} mode.")
        
        try:
            self._create_install_dir()
            
            # Skip interactive prompts in headless mode
            if not self.headless:
                self._prompt_downloads()
                
            self._install_python_packages()
            self._install_system_packages()
            
            # Use environment variables or command line args in headless mode
            self._setup_plastz()
            self._setup_blast_db()
            self._update_paths_in_script()
            self._verify_installations()
            
            logger.info("✅ Installation complete! ddPrimer is now ready to use.")
        except Exception as e:
            logger.error(f"❌ Setup failed: {str(e)}")
            sys.exit(1)
    
    def _create_install_dir(self):
        """Create installation directory if it doesn't exist."""
        if not os.path.exists(self.install_dir):
            os.makedirs(self.install_dir)
            logger.info(f"Created installation directory at {self.install_dir}")
        
        # Create support directory if on macOS
        if self.os_type == "Darwin" and not os.path.exists(SUPPORT_DIR):
            try:
                os.makedirs(SUPPORT_DIR, exist_ok=True)
                logger.info(f"Created support directory at {SUPPORT_DIR}")
            except PermissionError:
                logger.warning(f"Could not create {SUPPORT_DIR} due to permission error.")
                logger.warning("Some paths may need to be adjusted manually in ddPrimer.py")
    
    def _prompt_downloads(self):
        """Prompt user to manually download necessary files before proceeding."""
        logger.info("🔗 Please download the following before proceeding:")
        logger.info("1. NUPACK: https://nupack.org/download/")
        logger.info("2. MultiZ: https://www.bx.psu.edu/miller_lab/")
        logger.info("3. A FASTA file to create a BLAST+ database (if not already set up).")
        input("Press Enter once you have downloaded these files...")
    
    def _select_folder(self, prompt):
        """Open a file dialog to select a folder, or use command-line input if no GUI.
        
        Args:
            prompt (str): Dialog title prompt
            
        Returns:
            str: Selected folder path
        """
        # Check for default path in headless mode
        if self.headless:
            # Get the environment variable that corresponds to the prompt
            env_var = None
            if "MultiZ" in prompt:
                env_var = os.environ.get('DDPRIMER_MULTIZ_PATH')
            elif "PLastZ" in prompt:
                env_var = os.environ.get('DDPRIMER_PLASTZ_DIR')
                
            # If no environment variable, use a default in the install directory
            if not env_var:
                default_dir = os.path.join(self.install_dir, prompt.split()[1].lower())
                logger.info(f"Running in headless mode. Using default path: {default_dir}")
                return default_dir
            else:
                logger.info(f"Using path from environment variable: {env_var}")
                return env_var
        
        # Try to use GUI if available
        try:
            if not self.headless:
                import tkinter as tk
                from tkinter import filedialog
                root = tk.Tk()
                root.withdraw()
                folder_path = filedialog.askdirectory(title=prompt)
                if folder_path:
                    return folder_path
        except (ImportError, tk.TclError, Exception) as e:
            logger.info(f"GUI dialog not available ({str(e)}). Using command line input.")
        
        # Command line fallback
        print(f"\n{prompt}")
        folder_path = input("Enter the full path to the folder: ").strip()
        
        if not folder_path or not os.path.isdir(folder_path):
            logger.error(f"❌ Invalid directory path: {folder_path}")
            logger.info("Please enter a valid directory path.")
            return self._select_folder(prompt)  # Recursively prompt again
            
        return folder_path
    
    def _select_file(self, prompt, filetypes):
        """Open a file dialog to select a file, or use command-line input if no GUI.
        
        Args:
            prompt (str): Dialog title prompt
            filetypes (list): List of file type tuples (description, pattern)
            
        Returns:
            str: Selected file path
        """
        # Check for default path in headless mode
        if self.headless:
            # Get the environment variable that corresponds to the prompt
            env_var = None
            if "PLastZ" in prompt:
                env_var = os.environ.get('DDPRIMER_PLASTZ_PATH')
            elif "FASTA" in prompt:
                env_var = os.environ.get('DDPRIMER_FASTA_PATH')
                
            # If no environment variable, use a default in the install directory
            if not env_var:
                # Create a placeholder file if we're in headless mode
                if "PLastZ" in prompt:
                    default_file = os.path.join(self.install_dir, PLASTZ_DIR, PLASTZ_SCRIPT)
                    os.makedirs(os.path.dirname(default_file), exist_ok=True)
                    # Create an empty file if it doesn't exist
                    if not os.path.exists(default_file):
                        with open(default_file, 'w') as f:
                            f.write("# Placeholder PLastZ script - replace with actual script\n")
                elif "FASTA" in prompt:
                    default_file = os.path.join(self.install_dir, "Blast DB", "placeholder.fasta")
                    os.makedirs(os.path.dirname(default_file), exist_ok=True)
                    # Create a minimal FASTA file if it doesn't exist
                    if not os.path.exists(default_file):
                        with open(default_file, 'w') as f:
                            f.write(">placeholder\nACGT\n")
                else:
                    default_file = os.path.join(self.install_dir, "placeholder.txt")
                
                logger.info(f"Running in headless mode. Using default file: {default_file}")
                return default_file
            else:
                logger.info(f"Using file from environment variable: {env_var}")
                return env_var
        
        # Try to use GUI if available
        try:
            if not self.headless:
                import tkinter as tk
                from tkinter import filedialog
                root = tk.Tk()
                root.withdraw()
                file_path = filedialog.askopenfilename(title=prompt, filetypes=filetypes)
                if file_path:
                    return file_path
        except (ImportError, tk.TclError, Exception) as e:
            logger.info(f"GUI dialog not available ({str(e)}). Using command line input.")
        
        # Format filetypes for display
        extensions = []
        for desc, exts in filetypes:
            extensions.extend([ext.strip('*') for ext in exts.split()])
        
        # Command line fallback
        print(f"\n{prompt}")
        print(f"Acceptable file types: {', '.join(extensions)}")
        file_path = input("Enter the full path to the file: ").strip()
        
        if not file_path or not os.path.isfile(file_path):
            logger.error(f"❌ Invalid file path: {file_path}")
            logger.info("Please enter a valid file path.")
            return self._select_file(prompt, filetypes)  # Recursively prompt again
        
        # Basic extension check
        valid_extension = False
        for ext in extensions:
            if file_path.lower().endswith(ext.lower()):
                valid_extension = True
                break
        
        if not valid_extension:
            logger.warning(f"⚠️ File doesn't have an expected extension. Acceptable types: {', '.join(extensions)}")
            proceed = input("Continue anyway? (y/n): ").lower()
            if proceed != 'y':
                return self._select_file(prompt, filetypes)  # Recursively prompt again
            
        return file_path
    
    def _install_python_packages(self):
        """Install required Python packages using pip."""
        logger.info("📦 Installing required Python packages...")
        for package in PYTHON_PACKAGES:
            try:
                __import__(package)
                logger.info(f"✓ {package} is already installed.")
            except ImportError:
                logger.info(f"Installing {package}...")
                subprocess.run(
                    [sys.executable, "-m", "pip", "install", package], 
                    check=True, 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE
                )
                logger.info(f"✓ {package} installed successfully.")
    
    def _install_system_packages(self):
        """Install required system dependencies based on the operating system."""
        logger.info(f"🔧 Installing dependencies for {self.os_type}...")
        
        if self.os_type == "Linux":
            self._install_linux_packages()
            self._compile_multiz()
        elif self.os_type == "Darwin":  # macOS
            self._install_macos_packages()
            self._compile_multiz()
        elif self.os_type == "Windows":
            self._install_windows_packages()
        else:
            logger.error("❌ Unsupported operating system.")
            sys.exit(1)
    
    def _install_linux_packages(self):
        """Install system packages on Linux."""
        try:
            # Check if we have sudo access
            if self.headless:
                # In headless mode on HPC, we likely don't have sudo access
                logger.warning("Running in headless mode. Skipping system package installation.")
                logger.warning("Please ensure the following packages are installed:")
                for package in SYSTEM_PACKAGES["Linux"].values():
                    logger.warning(f"- {package}")
                return
                
            # Try to use apt-get
            subprocess.run(
                ["sudo", "apt-get", "update"], 
                check=True,
                stdout=subprocess.PIPE
            )
            
            packages = list(SYSTEM_PACKAGES["Linux"].values())
            logger.info(f"Installing packages: {', '.join(packages)}")
            
            subprocess.run(
                ["sudo", "apt-get", "install", "-y"] + packages, 
                check=True,
                stdout=subprocess.PIPE
            )
            logger.info("✔ System packages installed successfully.")
        except subprocess.SubprocessError as e:
            logger.error(f"❌ Failed to install system packages: {e}")
            logger.error("Please install them manually.")
    
    def _install_macos_packages(self):
        """Install system packages on macOS using Homebrew."""
        try:
            # Check if Homebrew is installed
            subprocess.run(
                ["brew", "--version"], 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                check=True
            )
        except (subprocess.SubprocessError, FileNotFoundError):
            logger.error("⚠ Homebrew not found! Please install it from https://brew.sh/")
            return
        
        # Install each required package
        for tool, package in SYSTEM_PACKAGES["Darwin"].items():
            if shutil.which(tool):
                logger.info(f"✔ {tool} is already installed.")
            else:
                logger.info(f"Installing {package} via Homebrew...")
                try:
                    subprocess.run(
                        ["brew", "install", package], 
                        check=True,
                        stdout=subprocess.PIPE
                    )
                    logger.info(f"✔ {package} installed successfully.")
                except subprocess.SubprocessError as e:
                    logger.error(f"❌ Failed to install {package}: {e}")
    
    def _install_windows_packages(self):
        """Provide instructions and help install packages on Windows."""
        logger.info("\n🔧 For Windows, please follow these steps:")
        logger.info("1. LastZ: Download from the official site and compile it using Cygwin/MSYS2.")
        logger.info("2. Samtools: Download from the official site and install it")
        logger.info("3. MultiZ: Follow the setup instructions for Cygwin/MSYS2.")
        logger.info("4. Primer3: Follow the manual compilation instructions from Primer3's official repository.")
        logger.info("5. BLAST+: Download the installer from the NCBI website and add it to your PATH.")
        
        if not self.headless:
            self._install_blast_windows()
    
    def _install_blast_windows(self):
        """Install BLAST+ on Windows using the provided installer."""
        logger.info("\n📦 Downloading and installing BLAST+ on Windows...")
        
        # Latest stable version URL
        blast_url = "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/blast-2.13.0+-win64.exe"
        blast_installer = Path("blast_installer.exe")
        
        try:
            # Download the BLAST+ installer
            logger.info(f"Downloading from {blast_url}...")
            urllib.request.urlretrieve(blast_url, blast_installer)
            logger.info(f"✔ Downloaded {blast_installer}.")
            
            # Run the installer
            logger.info("Running installer (silent mode)...")
            subprocess.run([str(blast_installer), "/S"], check=True)
            
            # Add BLAST+ to PATH
            blast_bin_path = r"C:\Program Files\NCBI\blast-2.13.0+\bin"
            os.environ["PATH"] += os.pathsep + blast_bin_path
            
            logger.info(f"✔ BLAST+ installed and added to PATH: {blast_bin_path}")
            logger.info("You may need to restart your terminal for the PATH changes to take effect.")
        except Exception as e:
            logger.error(f"❌ Failed to install BLAST+: {e}")
            logger.error("Please install it manually from the NCBI website.")
    
    def _compile_multiz(self):
        """Prompt user to select MultiZ folder and compile it."""
        logger.info("\n⚙️ Compiling MultiZ...")
        
        if self.headless and os.environ.get('DDPRIMER_MULTIZ_PATH'):
            # In headless mode with environment variable, just use the specified path
            self.multiz_path = os.environ.get('DDPRIMER_MULTIZ_PATH')
            logger.info(f"Using MultiZ from environment variable: {self.multiz_path}")
            return
        
        multiz_folder = self._select_folder("Select the MultiZ folder (where you extracted it)")
        self.multiz_path = multiz_folder
        
        # Create the destination directory for MultiZ
        multiz_dest = os.path.join(SUPPORT_DIR, "Multiz") if self.os_type == "Darwin" else os.path.join(self.install_dir, "Multiz")
        os.makedirs(multiz_dest, exist_ok=True)
        
        # If we're in headless mode and the folder doesn't exist or is a placeholder, skip compilation
        if self.headless and (not os.path.exists(multiz_folder) or 
                              not os.path.isfile(os.path.join(multiz_folder, "Makefile"))):
            logger.warning("MultiZ source not available in headless mode. Skipping compilation.")
            return
            
        makefile_path = os.path.join(multiz_folder, "Makefile")
        
        try:
            # Modify the Makefile to remove `-Werror`
            with open(makefile_path, "r") as f:
                makefile_content = f.readlines()
            
            modified_content = [line.replace("-Werror", "") for line in makefile_content]
            
            with open(makefile_path, "w") as f:
                f.writelines(modified_content)
            
            logger.info("✅ Removed `-Werror` from Makefile.")
            
            # Compile MultiZ
            logger.info("🔨 Running `make` in MultiZ directory...")
            subprocess.run(
                ["make"], 
                cwd=multiz_folder, 
                check=True,
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE
            )
            
            # Copy the executable to our destination
            multiz_binary = os.path.join(multiz_folder, "multiz")
            multiz_dest_bin = os.path.join(multiz_dest, "multiz")
            
            if os.path.exists(multiz_binary):
                shutil.copy2(multiz_binary, multiz_dest_bin)
                os.chmod(multiz_dest_bin, 0o755)  # Make executable
                logger.info(f"✔ MultiZ compiled and installed to {multiz_dest_bin}")
            else:
                logger.error("❌ MultiZ binary not found after compilation.")
                
        except Exception as e:
            logger.error(f"❌ Failed to compile MultiZ: {e}")
            logger.error("Please compile it manually according to the documentation.")
    
    def _setup_plastz(self):
        """Extract and set up the PLastZ implementation."""
        logger.info("\n🔄 Setting up PLastZ...")
        
        # Create the destination directory
        plastz_dest = os.path.join(SUPPORT_DIR, PLASTZ_DIR) if self.os_type == "Darwin" else os.path.join(self.install_dir, PLASTZ_DIR)
        os.makedirs(plastz_dest, exist_ok=True)
        
        # If we're in headless mode with an environment variable for PLastZ
        if self.headless and os.environ.get('DDPRIMER_PLASTZ_PATH'):
            plastz_script = os.environ.get('DDPRIMER_PLASTZ_PATH')
            logger.info(f"Using PLastZ from environment variable: {plastz_script}")
        else:
            # Ask user to provide the PLastZ script or use default
            if not self.headless:
                logger.info("Please select the PLastZ Python script (PLastZ.py)")
                plastz_script = self._select_file(
                    "Select PLastZ.py file", 
                    [("Python files", "*.py"), ("All files", "*.*")]
                )
            else:
                # In headless mode, create a placeholder
                plastz_script = os.path.join(plastz_dest, PLASTZ_SCRIPT)
                with open(plastz_script, 'w') as f:
                    f.write("# Placeholder PLastZ script - replace with actual script\n")
                logger.warning("Created placeholder PLastZ script. Replace with actual script.")
        
        try:
            # Copy the script
            dest_path = os.path.join(plastz_dest, PLASTZ_SCRIPT)
            shutil.copy2(plastz_script, dest_path)
            os.chmod(dest_path, 0o755)  # Make executable
            logger.info(f"✔ Copied {plastz_script} to {dest_path}")
            self.plastz_path = dest_path
                
            logger.info("✔ PLastZ setup completed.")
            
        except Exception as e:
            logger.error(f"❌ Failed to set up PLastZ: {e}")
    
    def _setup_blast_db(self):
        """Prompt user to select a FASTA file and create a BLAST+ database."""
        logger.info("\n🔬 Setting up BLAST+ database...")
        
        # Check if makeblastdb is available
        if not shutil.which("makeblastdb"):
            logger.error("❌ makeblastdb command not found. Please install BLAST+ first.")
            return
        
        # If in headless mode, check if we have a FASTA path from environment
        if self.headless and os.environ.get('DDPRIMER_FASTA_PATH'):
            fasta_file = os.environ.get('DDPRIMER_FASTA_PATH')
            logger.info(f"Using FASTA from environment variable: {fasta_file}")
        else:
            if not self.headless:
                fasta_file = self._select_file(
                    "Select the FASTA file for the BLAST+ database", 
                    [("FASTA files", "*.fasta *.fa *.fna")]
                )
            else:
                # In headless mode, create a minimal placeholder
                fasta_file = os.path.join(self.install_dir, "Blast DB", "placeholder.fasta")
                os.makedirs(os.path.dirname(fasta_file), exist_ok=True)
                with open(fasta_file, 'w') as f:
                    f.write(">placeholder\nACGT\n")
                logger.warning("Created placeholder FASTA file. Replace with actual file.")
        
        # Create the destination directory
        db_dest = os.path.join(SUPPORT_DIR, "Tair DB") if self.os_type == "Darwin" else os.path.join(self.install_dir, "Blast DB")
        os.makedirs(db_dest, exist_ok=True)
        
        db_name = os.path.splitext(os.path.basename(fasta_file))[0]
        db_path = os.path.join(db_dest, f"{db_name}_db")
        self.blast_db_path = db_path
        
        logger.info(f"🛠 Creating BLAST+ database at {db_path}...")
        try:
            subprocess.run(
                ["makeblastdb", "-in", fasta_file, "-dbtype", "nucl", "-out", db_path], 
                check=True,
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE
            )
            
            # Verify BLAST+ database setup
            if os.path.exists(db_path + ".nsq"):
                logger.info("✔ BLAST+ database setup complete.")
            else:
                logger.error("❌ Error: BLAST+ database creation failed.")
                return
            
            logger.info("\n🔧 To use this database, set the BLASTDB environment variable:")
            logger.info(f"   export BLASTDB={os.path.dirname(db_path)}")
        except subprocess.SubprocessError as e:
            logger.error(f"❌ Failed to create BLAST+ database: {e}")
    
    def _update_paths_in_script(self):
        """Update paths in the ddPrimer.py script to match the installation."""
        logger.info("\n🔧 Updating paths in ddPrimer.py...")
        
        # Find the ddPrimer.py script in the current directory or parent directory
        script_path = None
        current_dir = os.path.dirname(os.path.abspath(__file__))
        
        for root, _, files in os.walk(current_dir):
            for file in files:
                if file == "ddPrimer.py":
                    script_path = os.path.join(root, file)
                    break
            if script_path:
                break
        
        if not script_path:
            logger.warning("❌ Could not find ddPrimer.py to update paths.")
            return
        
        try:
            with open(script_path, 'r') as f:
                content = f.read()
            
            # Replace paths with installed locations
            path_updates = {}
            
            # PLastZ path
            if self.plastz_path:
                path_updates["PLASTZ_PATH"] = f'"{self.plastz_path}"'
            
            # MultiZ path
            if self.multiz_path:
                multiz_bin = os.path.join(os.path.dirname(self.multiz_path), "Multiz", "multiz")
                if os.path.exists(multiz_bin):
                    path_updates["MULTIZ_PATH"] = f'"{multiz_bin}"'
            
            # BLAST paths
            blastn_path = shutil.which("blastn")
            if blastn_path:
                path_updates["BLASTN_PATH"] = f'"{blastn_path}"'
            
            # BLAST DB path
            if self.blast_db_path:
                path_updates["DB_PATH"] = f'"{self.blast_db_path}"'
            
            # Primer3 path
            primer3_path = shutil.which("primer3_core")
            if primer3_path:
                path_updates["PRIMER3_CORE_PATH"] = f'"{primer3_path}"'
            
            # Apply the updates
            for var, new_path in path_updates.items():
                pattern = f'{var} = .*'
                replacement = f'{var} = {new_path}'
                content = re.sub(pattern, replacement, content)
            
            # Write the updated content back
            with open(script_path, 'w') as f:
                f.write(content)
            
            logger.info("✔ Updated paths in ddPrimer.py")
            
        except Exception as e:
            logger.error(f"❌ Failed to update paths in ddPrimer.py: {e}")
            logger.error("You may need to manually update the paths in the script.")
    
    def _verify_installations(self):
        """Verify installations of all dependencies."""
        logger.info("\n🔍 Verifying installations...")
        
        # Check command-line tools
        tools = {
            "primer3_core": "Primer3", 
            "lastz": "LastZ", 
            "blastn": "BLAST+",
            "samtools": "Samtools"
        }
        
        for cmd, name in tools.items():
            path = shutil.which(cmd)
            if path:
                logger.info(f"✔ {name} found at: {path}")
            else:
                logger.warning(f"❌ {name} not found in PATH.")
        
        # Verify NUPACK installation
        try:
            import nupack
            logger.info("✔ NUPACK is installed.")
        except ImportError:
            logger.warning("❌ NUPACK is not installed.")
        
        # Check PLastZ
        if self.plastz_path and os.path.exists(self.plastz_path):
            logger.info(f"✔ PLastZ script found at: {self.plastz_path}")
        else:
            logger.warning("❌ PLastZ script not found.")
        
        # Check if BLAST database was created
        if self.blast_db_path and os.path.exists(self.blast_db_path + ".nsq"):
            logger.info(f"✔ BLAST+ database exists at: {self.blast_db_path}")
        else:
            logger.warning("❌ BLAST+ database not verified.")
        
        # Check the MultiZ installation
        multiz_bin = os.path.join(os.path.dirname(self.multiz_path), "Multiz", "multiz")
        if os.path.exists(multiz_bin):
            logger.info(f"✔ MultiZ binary found at: {multiz_bin}")
        else:
            logger.warning("❌ MultiZ binary not found.")
        
        logger.info("\n📝 Final Configuration Instructions:")
        logger.info("1. Verify all paths in ddPrimer.py match your system")
        logger.info("2. Ensure all dependencies are correctly installed")
        logger.info("3. If any paths are incorrect, update them manually in the script")
        logger.info("4. Run `python ddPrimer.py` to start the primer design pipeline")

if __name__ == "__main__":
    setup = DdPrimerSetup()
    setup.run()