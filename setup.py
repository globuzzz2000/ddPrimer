from setuptools import setup, find_packages

setup(
    name="ddprimer",
    version="0.1.0",
    description="A pipeline for designing primers optimized for droplet digital PCR with cross-species support",
    author="Your Name",
    author_email="your.email@example.com",
    url="https://github.com/yourusername/ddPrimer",
    packages=find_packages(),
    install_requires=[
        "biopython",
        "pandas",
        "numpy",
        "tqdm",
        "primer3-py"
    ],
    entry_points={
        'console_scripts': [
            'ddprimer=ddprimer.pipeline:run_pipeline',
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    python_requires=">=3.7",
)

def main():
    try:
        success = run_pipeline()
        if success:
            logger.info("\n=== Pipeline completed successfully ===")
        else:
            logger.info("\n=== Pipeline completed with no results ===")
        return 0
    except Exception as e:
        logger.error(f"\n!!! Pipeline failed: {str(e)}")
        traceback.print_exc()
        return 1
    finally:
        # Redirect stderr to suppress the C++ cleanup messages
        import sys, os
        sys.stderr = open(os.devnull, 'w')