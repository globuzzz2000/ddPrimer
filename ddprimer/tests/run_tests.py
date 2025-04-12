#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test runner for ddPrimer - runs tests with coverage reporting
"""

import os
import sys
import subprocess
import argparse


def parse_args():
    """Parse command-line arguments for test runner."""
    parser = argparse.ArgumentParser(description='Run tests for ddPrimer with coverage reporting')
    
    parser.add_argument('--unit', action='store_true', help='Run only unit tests')
    parser.add_argument('--integration', action='store_true', help='Run only integration tests')
    parser.add_argument('--all', action='store_true', help='Run all tests (default)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    parser.add_argument('--html', action='store_true', help='Generate HTML coverage report')
    parser.add_argument('--xml', action='store_true', help='Generate XML coverage report')
    
    args = parser.parse_args()
    
    # Default to all tests if no specific test type is selected
    if not (args.unit or args.integration):
        args.all = True
    
    return args


def run_tests(args):
    """Run the specified tests with coverage reporting."""
    # Base pytest command
    cmd = ['pytest']
    
    # Add verbosity if requested
    if args.verbose:
        cmd.append('-v')
    
    # Add coverage
    cmd.extend(['--cov=ddPrimer', '--cov-report=term'])
    
    # Add HTML report if requested
    if args.html:
        cmd.append('--cov-report=html')
    
    # Add XML report if requested
    if args.xml:
        cmd.append('--cov-report=xml')
    
    # Add test paths based on selected test types
    if args.all:
        cmd.append('tests/')
    else:
        if args.unit:
            cmd.extend(['tests/test_utils/', 'tests/test_core/'])
        if args.integration:
            cmd.extend(['tests/test_modes/', 'tests/test_pipeline.py', 'tests/test_common.py'])
    
    # Print the command being run
    print(f"Running: {' '.join(cmd)}")
    
    # Run the tests
    result = subprocess.run(cmd)
    
    # Return the exit code
    return result.returncode


def main():
    """Main entry point for the test runner."""
    args = parse_args()
    
    # Make sure we're in the project root directory
    project_root = os.path.dirname(os.path.abspath(__file__))
    os.chdir(project_root)
    
    # Run the tests
    exit_code = run_tests(args)
    
    # Exit with the test result code
    sys.exit(exit_code)


if __name__ == '__main__':
    main()