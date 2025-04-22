#!/usr/bin/env python3
"""
MSP Scraper Wrapper

This is a wrapper script for backwards compatibility with the original msp_scraper.py script.
It simply imports and calls the main function from mspeeps.cli.

Usage:
    python msp_scraper_wrapper.py input_file.tsv [--output_dir OUTPUT_DIR] [--log_file LOG_FILE]

"""

import sys
from mspeeps.cli import main

if __name__ == "__main__":
    sys.exit(main()) 