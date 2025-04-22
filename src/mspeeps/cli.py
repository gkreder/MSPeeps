"""
Command-line interface for MSPeeps.

This module provides the command-line functionality for the package.
"""

import os
import sys
import argparse
import logging
import pandas as pd
from typing import List, Optional

from .utils import parse_input_file, setup_logging
from .core import process_file

def parse_args(args: Optional[List[str]] = None) -> argparse.Namespace:
    """
    Parse command-line arguments.
    
    Args:
        args: Command-line arguments (if None, uses sys.argv)
    
    Returns:
        Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="MSPeeps - Extract and convert mass spectrometry spectra to MSP format"
    )
    
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to the input TSV/Excel file"
    )
    
    parser.add_argument(
        "--output_dir",
        type=str,
        default="output",
        help="Directory to store output MSP files (default: 'output')"
    )
    
    parser.add_argument(
        "--log_file",
        type=str,
        default=None,
        help="Path to store log file (default: None, log to console only)"
    )
    
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    return parser.parse_args(args)

def main(args: Optional[List[str]] = None) -> int:
    """
    Main entry point for the command-line interface.
    
    Args:
        args: Command-line arguments (if None, uses sys.argv)
    
    Returns:
        Exit code (0 for success, non-zero for failure)
    """
    # Parse arguments
    parsed_args = parse_args(args)
    
    # Set up logging
    log_level = logging.DEBUG if parsed_args.verbose else logging.INFO
    logger = setup_logging(parsed_args.log_file, log_level)
    
    try:
        # Parse input file
        logger.info(f"Parsing input file: {parsed_args.input_file}")
        df = parse_input_file(parsed_args.input_file)
        logger.info(f"Found {len(df)} rows in input file")
        
        # Create output directory if it doesn't exist
        if not os.path.exists(parsed_args.output_dir):
            os.makedirs(parsed_args.output_dir)
            logger.info(f"Created output directory: {parsed_args.output_dir}")
        
        # Process each row
        success_count = 0
        error_count = 0
        
        for i, row in df.iterrows():
            try:
                logger.info(f"Processing row {i+1}/{len(df)}: {row.get('Molecule_name', 'Unknown')}")
                output_path = process_file(row, parsed_args.output_dir)
                
                if output_path:
                    success_count += 1
                    logger.info(f"Successfully processed row {i+1}: {os.path.basename(output_path)}")
                else:
                    error_count += 1
                    logger.error(f"Failed to process row {i+1}")
            except Exception as e:
                error_count += 1
                logger.error(f"Error processing row {i+1}: {str(e)}")
                if parsed_args.verbose:
                    logger.exception(e)
        
        # Summary
        logger.info(f"Processing complete: {success_count} succeeded, {error_count} failed")
        
        return 0 if error_count == 0 else 1
    
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        if parsed_args.verbose:
            logger.exception(e)
        return 1

if __name__ == "__main__":
    sys.exit(main()) 