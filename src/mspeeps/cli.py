"""
Command-line interface for MSPeeps.

This module provides the command-line functionality for the package.
"""

import os
import sys
import argparse
import logging
import json
import numpy as np
import pandas as pd
from typing import List, Optional

from .utils import parse_input_file, setup_logging, convert_smiles_to_inchi
from .core import process_file, extract_spectrum, format_msp
from .formula_matching import match_formula

def create_parser() -> argparse.ArgumentParser:
    """
    Create the command-line parser with all subcommands.
    
    Returns:
        Configured argument parser
    """
    parser = argparse.ArgumentParser(
        description="MSPeeps - Extract and convert mass spectrometry spectra to MSP format"
    )
    
    # Add version argument
    parser.add_argument(
        "--version",
        action="store_true",
        help="Show version and exit"
    )
    
    # Create subparsers for each command
    subparsers = parser.add_subparsers(dest="command", help="Command to run")
    
    # Batch command (existing functionality)
    add_batch_parser(subparsers)
    
    # Extract command (extract spectrum from mzML)
    add_extract_parser(subparsers)
    
    # Convert command (convert to MSP format)
    add_convert_parser(subparsers)
    
    # Match formula command
    add_match_formula_parser(subparsers)
    
    # Convert SMILES command
    add_convert_smiles_parser(subparsers)
    
    # Info command
    add_info_parser(subparsers)
    
    return parser

def add_batch_parser(subparsers) -> None:
    """Add the batch command parser (original functionality)."""
    parser = subparsers.add_parser(
        "batch",
        help="Process batch of spectra from an input file"
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

def add_extract_parser(subparsers) -> None:
    """Add the extract command parser."""
    parser = subparsers.add_parser(
        "extract",
        help="Extract spectrum from mzML file"
    )
    
    parser.add_argument(
        "--mzml_file",
        type=str,
        required=True,
        help="Path to the mzML file"
    )
    
    spectrum_group = parser.add_mutually_exclusive_group(required=True)
    spectrum_group.add_argument(
        "--spectrum_index",
        type=int,
        help="Index of the spectrum to extract"
    )
    spectrum_group.add_argument(
        "--retention_time",
        type=float,
        help="Retention time in seconds to find closest spectrum"
    )
    
    parser.add_argument(
        "--ms_level",
        type=int,
        default=2,
        help="MS level to filter spectra (default: 2 for MS/MS)"
    )
    
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output file path (default: stdout as JSON)"
    )
    
    parser.add_argument(
        "--format",
        choices=["json", "csv", "tsv"],
        default="json",
        help="Output format (default: json)"
    )

def add_convert_parser(subparsers) -> None:
    """Add the convert command parser."""
    parser = subparsers.add_parser(
        "convert",
        help="Convert spectrum data to MSP format"
    )
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--json_file",
        type=str,
        help="Input JSON file with spectrum data"
    )
    input_group.add_argument(
        "--mz_file",
        type=str,
        help="File containing m/z values (one per line)"
    )
    
    parser.add_argument(
        "--intensity_file",
        type=str,
        help="File containing intensity values (required if --mz_file is used)"
    )
    
    parser.add_argument(
        "--metadata_file",
        type=str,
        help="JSON file with metadata (optional if --mz_file is used)"
    )
    
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output file path (default: stdout)"
    )
    
    parser.add_argument(
        "--intensity_cutoff",
        type=float,
        default=0,
        help="Raw intensity cutoff (default: 0)"
    )

def add_match_formula_parser(subparsers) -> None:
    """Add the match-formula command parser."""
    parser = subparsers.add_parser(
        "match-formula",
        help="Match m/z values to molecular formulas"
    )
    
    parser.add_argument(
        "--mz_values",
        type=str,
        required=True,
        help="Comma-separated m/z values or file with one value per line"
    )
    
    parser.add_argument(
        "--parent_formula",
        type=str,
        required=True,
        help="Molecular formula of the parent molecule"
    )
    
    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.002,
        help="Mass tolerance in Da (default: 0.002)"
    )
    
    parser.add_argument(
        "--charge",
        type=int,
        default=1,
        help="Charge state (default: 1)"
    )
    
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output file path (default: stdout)"
    )
    
    parser.add_argument(
        "--format",
        choices=["tsv", "csv", "json", "msp"],
        default="tsv",
        help="Output format (default: tsv)"
    )

def add_convert_smiles_parser(subparsers) -> None:
    """Add the convert-smiles command parser."""
    parser = subparsers.add_parser(
        "convert-smiles",
        help="Convert SMILES to InChI and InChIKey"
    )
    
    parser.add_argument(
        "--smiles",
        type=str,
        required=True,
        help="SMILES string to convert"
    )
    
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output file path (default: stdout)"
    )

def add_info_parser(subparsers) -> None:
    """Add the info command parser."""
    parser = subparsers.add_parser(
        "info",
        help="Display information about mzML files"
    )
    
    parser.add_argument(
        "--mzml_file",
        type=str,
        required=True,
        help="Path to the mzML file"
    )
    
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output file path (default: stdout)"
    )
    
    parser.add_argument(
        "--format",
        choices=["text", "json"],
        default="text",
        help="Output format (default: text)"
    )

def batch_command(args: argparse.Namespace, logger: logging.Logger) -> int:
    """
    Run the batch command (original functionality).
    
    Args:
        args: Command-line arguments
        logger: Logger instance
    
    Returns:
        Exit code
    """
    try:
        # Parse input file
        logger.info(f"Parsing input file: {args.input_file}")
        df = parse_input_file(args.input_file)
        logger.info(f"Found {len(df)} rows in input file")
        
        # Create output directory if it doesn't exist
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)
            logger.info(f"Created output directory: {args.output_dir}")
        
        # Process each row
        success_count = 0
        error_count = 0
        
        for i, row in df.iterrows():
            try:
                logger.info(f"Processing row {i+1}/{len(df)}: {row.get('Molecule_name', 'Unknown')}")
                output_path = process_file(row, args.output_dir)
                
                if output_path:
                    success_count += 1
                    logger.info(f"Successfully processed row {i+1}: {os.path.basename(output_path)}")
                else:
                    error_count += 1
                    logger.error(f"Failed to process row {i+1}")
            except Exception as e:
                error_count += 1
                logger.error(f"Error processing row {i+1}: {str(e)}")
                if args.verbose:
                    logger.exception(e)
        
        # Summary
        logger.info(f"Processing complete: {success_count} succeeded, {error_count} failed")
        
        return 0 if error_count == 0 else 1
    
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        if args.verbose:
            logger.exception(e)
        return 1

def extract_command(args: argparse.Namespace, logger: logging.Logger) -> int:
    """
    Run the extract command.
    
    Args:
        args: Command-line arguments
        logger: Logger instance
    
    Returns:
        Exit code
    """
    try:
        # Extract spectrum
        logger.info(f"Extracting spectrum from: {args.mzml_file}")
        mz_array, intensity_array, metadata = extract_spectrum(
            mzml_path=args.mzml_file,
            spectrum_index=args.spectrum_index,
            retention_time=args.retention_time,
            ms_level=args.ms_level
        )
        
        # Format output data
        output_data = {
            "mz_array": mz_array.tolist(),
            "intensity_array": intensity_array.tolist(),
            "metadata": metadata
        }
        
        # Output the data
        if args.format == "json":
            output_str = json.dumps(output_data, indent=2)
        elif args.format == "csv":
            output_str = "m/z,intensity\n"
            output_str += "\n".join(f"{mz},{intensity}" for mz, intensity in zip(mz_array, intensity_array))
        elif args.format == "tsv":
            output_str = "m/z\tintensity\n"
            output_str += "\n".join(f"{mz}\t{intensity}" for mz, intensity in zip(mz_array, intensity_array))
        
        # Write output
        if args.output:
            with open(args.output, "w") as f:
                f.write(output_str)
            logger.info(f"Output written to: {args.output}")
        else:
            print(output_str)
        
        return 0
    
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        logger.exception(e)
        return 1

def convert_command(args: argparse.Namespace, logger: logging.Logger) -> int:
    """
    Run the convert command.
    
    Args:
        args: Command-line arguments
        logger: Logger instance
    
    Returns:
        Exit code
    """
    try:
        # Get input data
        if args.json_file:
            # Load data from JSON file
            with open(args.json_file, "r") as f:
                data = json.load(f)
            
            mz_array = np.array(data["mz_array"])
            intensity_array = np.array(data["intensity_array"])
            metadata = data.get("metadata", {})
            
            # Use row_data from JSON if available, otherwise create empty
            row_data = pd.Series(data.get("row_data", {}))
            
        else:
            # Load m/z values from file
            with open(args.mz_file, "r") as f:
                mz_array = np.array([float(line.strip()) for line in f if line.strip()])
            
            # Load intensity values from file
            if not args.intensity_file:
                logger.error("--intensity_file is required when using --mz_file")
                return 1
                
            with open(args.intensity_file, "r") as f:
                intensity_array = np.array([float(line.strip()) for line in f if line.strip()])
            
            # Check array lengths match
            if len(mz_array) != len(intensity_array):
                logger.error(f"Mismatch in array lengths: {len(mz_array)} m/z values and {len(intensity_array)} intensity values")
                return 1
            
            # Load metadata if provided
            if args.metadata_file:
                with open(args.metadata_file, "r") as f:
                    metadata = json.load(f)
            else:
                metadata = {}
            
            # Create empty row_data
            row_data = pd.Series({})
        
        # Apply intensity cutoff
        if args.intensity_cutoff > 0:
            mask = intensity_array >= args.intensity_cutoff
            mz_array = mz_array[mask]
            intensity_array = intensity_array[mask]
            logger.info(f"Applied intensity cutoff {args.intensity_cutoff}: {len(mz_array)} peaks remain")
        
        # Format as MSP
        msp_data = format_msp(
            mz_array,
            intensity_array,
            metadata,
            row_data,
            raw_intensity_cutoff=0  # Already applied cutoff above
        )
        
        # Write output
        if args.output:
            with open(args.output, "w") as f:
                f.write(msp_data)
            logger.info(f"MSP data written to: {args.output}")
        else:
            print(msp_data)
        
        return 0
    
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        logger.exception(e)
        return 1

def match_formula_command(args: argparse.Namespace, logger: logging.Logger) -> int:
    """
    Run the match-formula command.
    
    Args:
        args: Command-line arguments
        logger: Logger instance
    
    Returns:
        Exit code
    """
    try:
        # Get m/z values
        if os.path.exists(args.mz_values):
            # Read from file
            with open(args.mz_values, "r") as f:
                mz_values = np.array([float(line.strip()) for line in f if line.strip()])
        else:
            # Parse comma-separated values
            mz_values = np.array([float(x.strip()) for x in args.mz_values.split(",") if x.strip()])
        
        logger.info(f"Matching {len(mz_values)} m/z values to formulas")
        
        # Match formulas
        matches = match_formula(
            mz_array=mz_values,
            parent_formula=args.parent_formula,
            tolerance_da=args.tolerance,
            charge=args.charge
        )
        
        # Format output
        if args.format == "tsv":
            output_str = "m/z\tformula\texact_mass\tmass_error\n"
            for mz, (formula, exact_mass, error) in zip(mz_values, matches):
                if formula and exact_mass:
                    output_str += f"{mz:.6f}\t{formula}\t{exact_mass:.6f}\t{error:.6f}\n"
                else:
                    output_str += f"{mz:.6f}\t\t\t\n"
        
        elif args.format == "csv":
            output_str = "m/z,formula,exact_mass,mass_error\n"
            for mz, (formula, exact_mass, error) in zip(mz_values, matches):
                if formula and exact_mass:
                    output_str += f"{mz:.6f},{formula},{exact_mass:.6f},{error:.6f}\n"
                else:
                    output_str += f"{mz:.6f},,,\n"
        
        elif args.format == "json":
            results = []
            for mz, (formula, exact_mass, error) in zip(mz_values, matches):
                if formula and exact_mass:
                    results.append({
                        "m/z": mz,
                        "formula": formula,
                        "exact_mass": exact_mass,
                        "mass_error": error
                    })
                else:
                    results.append({
                        "m/z": mz,
                        "formula": None,
                        "exact_mass": None,
                        "mass_error": None
                    })
            output_str = json.dumps(results, indent=2)
        
        elif args.format == "msp":
            output_str = "NUM PEAKS: {}\n".format(len(mz_values))
            for mz, (formula, exact_mass, error) in zip(mz_values, matches):
                if formula and exact_mass:
                    output_str += f"{mz:.6f} 1000 \"{formula}\" {exact_mass:.6f} {error:.6f}\n"
                else:
                    output_str += f"{mz:.6f} 1000\n"
        
        # Write output
        if args.output:
            with open(args.output, "w") as f:
                f.write(output_str)
            logger.info(f"Output written to: {args.output}")
        else:
            print(output_str)
        
        return 0
    
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        logger.exception(e)
        return 1

def convert_smiles_command(args: argparse.Namespace, logger: logging.Logger) -> int:
    """
    Run the convert-smiles command.
    
    Args:
        args: Command-line arguments
        logger: Logger instance
    
    Returns:
        Exit code
    """
    try:
        # Convert SMILES to InChI and InChIKey
        inchi, inchikey = convert_smiles_to_inchi(args.smiles)
        
        if not inchi:
            logger.error(f"Failed to convert SMILES: {args.smiles}")
            return 1
        
        # Format output
        output_data = {
            "SMILES": args.smiles,
            "InChI": inchi,
            "InChIKey": inchikey
        }
        
        output_str = json.dumps(output_data, indent=2)
        
        # Write output
        if args.output:
            with open(args.output, "w") as f:
                f.write(output_str)
            logger.info(f"Output written to: {args.output}")
        else:
            print(output_str)
        
        return 0
    
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        logger.exception(e)
        return 1

def info_command(args: argparse.Namespace, logger: logging.Logger) -> int:
    """
    Run the info command.
    
    Args:
        args: Command-line arguments
        logger: Logger instance
    
    Returns:
        Exit code
    """
    try:
        import pymzml
        
        # Check if file exists
        if not os.path.exists(args.mzml_file):
            logger.error(f"mzML file not found: {args.mzml_file}")
            return 1
        
        # Open mzML file
        logger.info(f"Analyzing mzML file: {args.mzml_file}")
        run = pymzml.run.Reader(args.mzml_file)
        
        # Get file info
        spectrum_count = 0
        ms1_count = 0
        ms2_count = 0
        ms_levels = set()
        rt_min = float('inf')
        rt_max = 0
        
        # Analyze first few spectra
        max_preview = 100
        preview_count = 0
        
        for spectrum in run:
            spectrum_count += 1
            
            # Get MS level
            ms_level = spectrum.ms_level
            ms_levels.add(ms_level)
            
            if ms_level == 1:
                ms1_count += 1
            elif ms_level == 2:
                ms2_count += 1
            
            # Get retention time
            if hasattr(spectrum, 'scan_time'):
                rt = spectrum.scan_time[0]
                rt_min = min(rt_min, rt)
                rt_max = max(rt_max, rt)
            
            # Limit full scan to max_preview spectra
            preview_count += 1
            if preview_count >= max_preview:
                break
        
        # Prepare info data
        info = {
            "file_path": args.mzml_file,
            "file_size_mb": os.path.getsize(args.mzml_file) / (1024 * 1024),
            "spectrum_count": spectrum_count,
            "ms1_count": ms1_count,
            "ms2_count": ms2_count,
            "ms_levels": sorted(list(ms_levels)),
            "retention_time_range": [rt_min, rt_max] if rt_min != float('inf') else None,
            "preview_limit": max_preview
        }
        
        # Format output
        if args.format == "json":
            output_str = json.dumps(info, indent=2)
        else:  # text format
            output_str = f"File: {info['file_path']}\n"
            output_str += f"Size: {info['file_size_mb']:.2f} MB\n"
            output_str += f"Spectra: {info['spectrum_count']} (preview of {info['preview_limit']})\n"
            output_str += f"MS1 Spectra: {info['ms1_count']}\n"
            output_str += f"MS2 Spectra: {info['ms2_count']}\n"
            output_str += f"MS Levels: {', '.join(str(x) for x in info['ms_levels'])}\n"
            
            if info['retention_time_range']:
                rt_min_min = info['retention_time_range'][0] / 60
                rt_max_min = info['retention_time_range'][1] / 60
                output_str += f"RT Range: {rt_min_min:.2f} - {rt_max_min:.2f} min ({info['retention_time_range'][0]:.2f} - {info['retention_time_range'][1]:.2f} sec)\n"
            else:
                output_str += "RT Range: Not available\n"
        
        # Write output
        if args.output:
            with open(args.output, "w") as f:
                f.write(output_str)
            logger.info(f"File info written to: {args.output}")
        else:
            print(output_str)
        
        return 0
    
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        logger.exception(e)
        return 1

def main(args: Optional[List[str]] = None) -> int:
    """
    Main entry point for the command-line interface.
    
    Args:
        args: Command-line arguments (if None, uses sys.argv)
    
    Returns:
        Exit code (0 for success, non-zero for failure)
    """
    # Parse arguments
    parser = create_parser()
    parsed_args = parser.parse_args(args)
    
    # Show version and exit if requested
    if hasattr(parsed_args, 'version') and parsed_args.version:
        from . import __version__
        print(f"MSPeeps version {__version__}")
        return 0
    
    # Check if a command was specified
    if not hasattr(parsed_args, 'command') or not parsed_args.command:
        parser.print_help()
        return 1
    
    # Set up logging
    log_level = logging.DEBUG if hasattr(parsed_args, 'verbose') and parsed_args.verbose else logging.INFO
    log_file = parsed_args.log_file if hasattr(parsed_args, 'log_file') else None
    logger = setup_logging(log_file, log_level)
    
    # Execute the appropriate command
    if parsed_args.command == "batch":
        return batch_command(parsed_args, logger)
    elif parsed_args.command == "extract":
        return extract_command(parsed_args, logger)
    elif parsed_args.command == "convert":
        return convert_command(parsed_args, logger)
    elif parsed_args.command == "match-formula":
        return match_formula_command(parsed_args, logger)
    elif parsed_args.command == "convert-smiles":
        return convert_smiles_command(parsed_args, logger)
    elif parsed_args.command == "info":
        return info_command(parsed_args, logger)
    else:
        logger.error(f"Unknown command: {parsed_args.command}")
        return 1

if __name__ == "__main__":
    sys.exit(main()) 