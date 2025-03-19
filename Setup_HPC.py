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
import argparse

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(message)s'
)
logger = logging.getLogger('ddprimer_setup')

# Required Python Packages
PYTHON_PACKAGES = ["pandas", "tqdm"]

# System Packages mapping: command_name -> package_name
SYSTEM_PACKAGES = {
    "Linux": {
        "primer3_core": "primer3",
        "blastn": "ncbi-blast+",
        "lastz": "lastz"
    },
    "Darwin": {  # macOS
        "primer3_core": "primer3",
        "blastn": "blast",
        "lastz": "lastz"
    },
    "Windows": {
        # Windows packages are handled differently
    }
}

class DdPrimerSetup:
    """Main class for setting up ddPrimer dependencies."""
    
    def __init__(self):
        """Initialize setup environment variables."""
        self.os_type = platform.system()
        self.blast_db_path = None
        self.multiz_path = None
    
    def run(self):
        """Run the complete setup process."""
        logger.info("=== Setting up ddPrimer Dependencies ===")
        
        try:
            self._prompt_downloads()
            self._install_python_packages()
            self._install_nupack()
            self._install_system_packages()
            self._setup_blast_db()
            self._verify_installations()
            
            logger.info("✅ Installation complete! ddPrimer is now ready to use.")
            logger.info("⚠️ Reminder: Please ensure paths for MultiZ, Primer3, Blast+, and BLAST database are correctly set in the main script.")
        except Exception as e:
            logger.error(f"❌ Setup failed: {str(e)}")
            sys.exit(1)
    
    def _prompt_downloads(self):
        """Prompt user to manually download necessary files before proceeding."""
        logger.info("🔗 Please download the following before proceeding:")
        logger.info("1. NUPACK: https://nupack.org/download/")
        logger.info("2. MultiZ: https://www.bx.psu.edu/miller_lab/")
        logger.info("3. A FASTA file to create a BLAST+ database (if not already set up).")
        input("Press Enter once you have downloaded these files...")
    
    def _select_folder(self, prompt):
        folder_path = input(f"{prompt} (Enter full path): ").strip()
        if not os.path.isdir(folder_path):
            logger.error(f"❌ Directory not found: {folder_path}")
            sys.exit(1)
        return folder_path
    
    def _select_file(self, prompt):
        file_path = input(f"{prompt} (Enter full path): ").strip()
        if not os.path.isfile(file_path):
            logger.error(f"❌ File not found: {file_path}")
            sys.exit(1)
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
    
    def _install_nupack(self):
        """Prompt user to select the NUPACK .whl file and install it."""
        logger.info("📦 Installing NUPACK...")
        
        # First check if NUPACK is already installed
        try:
            import nupack
            logger.info("✓ NUPACK is already installed.")
            return
        except ImportError:
            pass
        
        whl_file = input("Enter the full path to the NUPACK .whl file: ").strip()
        
        logger.info(f"Installing {whl_file}...")
        try:
            subprocess.run(
                [sys.executable, "-m", "pip", "install", whl_file], 
                check=True,
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE
            )
            
            # Verify NUPACK installation
            import nupack
            logger.info("✓ NUPACK installed successfully.")
        except (subprocess.SubprocessError, ImportError) as e:
            logger.error(f"❌ Error: NUPACK installation failed: {e}")
            logger.error("Please install it manually.")
    
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
        """Check for required tools in PATH instead of installing system packages."""
        logger.info("🔍 Checking for required system tools...")
        missing_tools = []
        
        for tool in SYSTEM_PACKAGES["Linux"]:
            if shutil.which(tool):
                logger.info(f"✔ {tool} is already installed.")
            else:
                logger.warning(f"❌ {tool} not found. Please install it manually.")
                missing_tools.append(tool)
        
        if missing_tools:
            logger.error(f"❌ The following tools are missing: {', '.join(missing_tools)}")
            logger.error("Please load them via 'module load' or install them manually.")
    
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
        logger.info("   - Cygwin/MSYS2 setup: Follow the instructions at https://www.cygwin.com/ or https://www.msys2.org/")
        logger.info("   - After installing, use the following commands to install LastZ:")
        logger.info("     `pacman -Syu`")
        logger.info("     `pacman -S lastz`")
        logger.info("2. MultiZ: Follow the setup instructions for Cygwin/MSYS2.")
        logger.info("3. Primer3: Follow the manual compilation instructions from Primer3's official repository.")
        logger.info("4. BLAST+: Download the installer from the NCBI website and add it to your PATH.")
        
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
        multiz_folder = input("Enter the full path to the MultiZ folder: ").strip()
        self.multiz_path = multiz_folder
        
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
            
            logger.info("✔ MultiZ compiled successfully.")
        except Exception as e:
            logger.error(f"❌ Failed to compile MultiZ: {e}")
            logger.error("Please compile it manually according to the documentation.")
    
    def _setup_blast_db(self):
        """Prompt user to select a FASTA file and create a BLAST+ database."""
        logger.info("\n🔬 Setting up BLAST+ database...")
        
        # Check if makeblastdb is available
        if not shutil.which("makeblastdb"):
            logger.error("❌ makeblastdb command not found. Please install BLAST+ first.")
            return
        
        fasta_file = input("Enter the full path to the FASTA file for BLAST+ database: ").strip()
        
        db_name = os.path.splitext(os.path.basename(fasta_file))[0] + "_db"
        db_path = os.path.join(os.path.dirname(fasta_file), db_name)
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
    
    def _verify_installations(self):
        """Verify installations of all dependencies."""
        logger.info("\n🔍 Verifying installations...")
        
        # Check command-line tools
        tools = ["primer3_core", "lastz", "blastn", "multiz"]
        for tool in tools:
            path = shutil.which(tool)
            if path:
                logger.info(f"✔ {tool} found at: {path}")
            else:
                logger.warning(f"❌ {tool} not found in PATH.")
        
        # Verify NUPACK installation
        try:
            import nupack
            logger.info("✔ NUPACK is installed.")
        except ImportError:
            logger.warning("❌ NUPACK is not installed.")
        
        # Check if BLAST database was created
        if self.blast_db_path and os.path.exists(self.blast_db_path + ".nsq"):
            logger.info(f"✔ BLAST+ database exists at: {self.blast_db_path}")
        else:
            logger.warning("❌ BLAST+ database not verified.")
        
        # Check the MultiZ installation
        if self.multiz_path:
            multiz_binary = os.path.join(self.multiz_path, "multiz")
            if os.path.exists(multiz_binary):
                logger.info(f"✔ MultiZ binary found at: {multiz_binary}")
            else:
                logger.warning("❌ MultiZ binary not found.")
        else:
            logger.warning("❌ MultiZ path not set.")


def main():
    """Main entry point for the setup script."""
    setup = DdPrimerSetup()
    setup.run()


if __name__ == "__main__":
    main()