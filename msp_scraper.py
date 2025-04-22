#!/usr/bin/env python3
"""
MSP File Scraper

A tool for scraping spectra from mzML files and converting them to MSP format.
"""

import os
import sys
import pandas as pd
import numpy as np
import logging
import argparse
from typing import Dict, List, Tuple, Optional, Union
import traceback
# Import fragfit for formula matching
import fragfit

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_input_file(file_path: str) -> pd.DataFrame:
    """Parse the input TSV/Excel file."""
    if file_path.endswith('.tsv'):
        df = pd.read_csv(file_path, sep='\t')
    elif file_path.endswith(('.xlsx', '.xls')):
        df = pd.read_excel(file_path)
    else:
        raise ValueError("Input file must be a TSV or Excel file.")
    
    # Check required columns
    required_cols = ['Molecule_name']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {', '.join(missing_cols)}")
    
    return df

def extract_spectrum(mzml_path: str, spectrum_index: Optional[int] = None, 
                     retention_time: Optional[float] = None, 
                     ms_level: int = 2) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Extract spectrum from mzML file using index or RT.
    
    Args:
        mzml_path: Path to the mzML file
        spectrum_index: Index of the spectrum to extract
        retention_time: Retention time to find closest spectrum (in seconds)
        ms_level: MS level to filter spectra (default: 2 for MS/MS)
    
    Returns:
        Tuple of (mz_array, intensity_array, metadata)
    """
    import pymzml
    
    # Check if file exists
    if not os.path.exists(mzml_path):
        raise FileNotFoundError(f"mzML file not found: {mzml_path}")
    
    # Open mzML file
    run = pymzml.run.Reader(mzml_path)
    
    # If spectrum index is provided, extract by index
    if spectrum_index is not None:
        # Adjust for 0-based indexing if needed (spec says Agilent uses 1-based)
        adjusted_index = spectrum_index
        
        # Try to get the spectrum
        try:
            # pymzML uses different access methods depending on version
            spectrum = None
            if hasattr(run, 'get_by_id'):
                spectrum = run.get_by_id(adjusted_index)
            else:
                # Fallback to iterating (less efficient)
                for i, spec in enumerate(run):
                    if i + 1 == adjusted_index:  # +1 for 1-based indexing
                        spectrum = spec
                        break
            
            if spectrum is None:
                raise ValueError(f"Spectrum index {spectrum_index} not found in file.")
            
            # Check MS level if provided
            if ms_level and spectrum.ms_level != ms_level:
                logger.warning(f"Spectrum has MS level {spectrum.ms_level}, expected {ms_level}")
            
            # Get m/z and intensity arrays
            mz_array = spectrum.mz
            intensity_array = spectrum.i
            
            # Get metadata
            metadata = {
                'ms_level': spectrum.ms_level,
                'retention_time': spectrum.scan_time[0] if hasattr(spectrum, 'scan_time') else None,
                'precursor_mz': spectrum.selected_precursors[0]['mz'] if hasattr(spectrum, 'selected_precursors') and spectrum.selected_precursors else None,
                'precursor_charge': spectrum.selected_precursors[0].get('charge', None) if hasattr(spectrum, 'selected_precursors') and spectrum.selected_precursors else None,
            }
            
            return mz_array, intensity_array, metadata
            
        except Exception as e:
            raise ValueError(f"Error extracting spectrum {spectrum_index}: {str(e)}")
    
    # If retention time is provided, find closest spectrum
    elif retention_time is not None:
        closest_spectrum = None
        min_rt_diff = float('inf')
        
        for spectrum in run:
            # Skip spectra with wrong MS level
            if ms_level and spectrum.ms_level != ms_level:
                continue
                
            # Get spectrum retention time
            if hasattr(spectrum, 'scan_time'):
                rt = spectrum.scan_time[0]
                rt_diff = abs(rt - retention_time)
                
                if rt_diff < min_rt_diff:
                    min_rt_diff = rt_diff
                    closest_spectrum = spectrum
        
        if closest_spectrum is None:
            raise ValueError(f"No spectrum found with MS level {ms_level} near retention time {retention_time}")
        
        # Get m/z and intensity arrays
        mz_array = closest_spectrum.mz
        intensity_array = closest_spectrum.i
        
        # Get metadata
        metadata = {
            'ms_level': closest_spectrum.ms_level,
            'retention_time': closest_spectrum.scan_time[0] if hasattr(closest_spectrum, 'scan_time') else None,
            'precursor_mz': closest_spectrum.selected_precursors[0]['mz'] if hasattr(closest_spectrum, 'selected_precursors') and closest_spectrum.selected_precursors else None,
            'precursor_charge': closest_spectrum.selected_precursors[0].get('charge', None) if hasattr(closest_spectrum, 'selected_precursors') and closest_spectrum.selected_precursors else None,
        }
        
        return mz_array, intensity_array, metadata
    
    else:
        raise ValueError("Either spectrum_index or retention_time must be provided")

def convert_smiles_to_inchi(smiles: str) -> Tuple[str, str]:
    """
    Convert SMILES to InChI and InChIKey.
    
    Args:
        smiles: SMILES string
    
    Returns:
        Tuple of (InChI, InChIKey)
    """
    if not smiles or pd.isna(smiles):
        return "", ""
    
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Could not parse SMILES: {smiles}")
            return "", ""
        
        inchi = Chem.MolToInchi(mol)
        inchikey = Chem.MolToInchiKey(mol)
        return inchi, inchikey
    except Exception as e:
        logger.warning(f"Error converting SMILES to InChI: {str(e)}")
        return "", ""

def match_formula(mz_array: np.ndarray, parent_formula: str, 
                 tolerance_da: float, charge: int = 1) -> List[Tuple[str, float, float]]:
    """
    Match peaks to possible fragment formulas.
    
    Args:
        mz_array: Array of m/z values
        parent_formula: Molecular formula of the parent molecule
        tolerance_da: Mass tolerance in Da
        charge: Charge state of the fragments (default: 1)
    
    Returns:
        List of tuples (formula, exact_mass, mass_error)
    """
    if not parent_formula or pd.isna(parent_formula):
        logger.warning("No parent formula provided for formula matching")
        return [(None, None, None)] * len(mz_array)
    
    try:
        results = []
        
        for mz in mz_array:
            # Use fragfit.find_best_form to find best matching formula
            formula, exact_mass, error = fragfit.find_best_form(
                mass=mz,
                parent_form=parent_formula,
                tolerance_da=tolerance_da,
                charge=charge
            )
            
            results.append((formula, exact_mass, error))
        
        return results
    except Exception as e:
        logger.error(f"Error in formula matching: {str(e)}")
        return [(None, None, None)] * len(mz_array)

def format_msp(mz_array: np.ndarray, intensity_array: np.ndarray, 
               metadata: Dict, row_data: pd.Series, 
               raw_intensity_cutoff: float = 0) -> str:
    """
    Format the data in MSP format.
    
    Args:
        mz_array: Array of m/z values
        intensity_array: Array of intensity values
        metadata: Additional metadata from the spectrum
        row_data: Row data from the input file
        raw_intensity_cutoff: Intensity cutoff for peaks
    
    Returns:
        Formatted MSP string
    """
    # Apply intensity cutoff
    if raw_intensity_cutoff > 0:
        mask = intensity_array >= raw_intensity_cutoff
        mz_array = mz_array[mask]
        intensity_array = intensity_array[mask]
    
    # Check if formula matching should be applied
    do_formula_matching = (pd.notna(row_data.get('Molecular_formula')) and 
                          pd.notna(row_data.get('Formula_Matching_Tolerance')))
    
    # Perform formula matching if needed
    formula_matches = None
    if do_formula_matching:
        try:
            # Get charge from metadata if available
            charge = metadata.get('precursor_charge', 1)
            if charge is None:
                charge = 1
            
            logger.info(f"Performing formula matching with parent formula: {row_data['Molecular_formula']}")
            formula_matches = match_formula(
                mz_array, 
                row_data['Molecular_formula'],
                float(row_data['Formula_Matching_Tolerance']),
                charge
            )
        except Exception as e:
            logger.warning(f"Error in formula matching: {str(e)}")
            do_formula_matching = False
    
    # Calculate number of peaks and total intensity
    num_peaks = len(mz_array)
    total_intensity = sum(intensity_array)
    
    # Start building MSP string
    msp_lines = []
    
    # Add NAME as the first field (required)
    msp_lines.append(f"NAME: {row_data['Molecule_name']}")
    
    # Add metadata fields from the input row
    for col, value in row_data.items():
        if pd.notna(value) and col != 'Molecule_name' and col != 'InChI' and col != 'InChIKey':
            # Skip the NAME field as we already added it
            # Convert column names to match MSP format (e.g., Molecular_formula -> FORMULA)
            field_name = col.upper().replace('_', ' ')
            msp_lines.append(f"{field_name}: {value}")
    
    # Add calculated InChI and InChIKey if SMILES is provided
    if pd.notna(row_data.get('SMILES')):
        inchi, inchikey = convert_smiles_to_inchi(row_data['SMILES'])
        if inchi:
            msp_lines.append(f"INCHI: {inchi}")
        if inchikey:
            msp_lines.append(f"INCHIKEY: {inchikey}")
    
    # Add spectrum metadata
    if metadata.get('retention_time') is not None:
        msp_lines.append(f"RETENTIONTIME: {metadata['retention_time']:.2f}")
    
    if metadata.get('precursor_mz') is not None:
        msp_lines.append(f"PRECURSORMZ: {metadata['precursor_mz']:.6f}")
    
    if metadata.get('ms_level') is not None:
        msp_lines.append(f"MSLEVEL: {metadata['ms_level']}")
    
    # Add spectrum info
    msp_lines.append(f"NUM PEAKS: {num_peaks}")
    
    # Add peak data
    peak_lines = []
    for i, (mz, intensity) in enumerate(zip(mz_array, intensity_array)):
        if do_formula_matching and formula_matches[i][0] is not None:
            formula, exact_mass, _ = formula_matches[i]
            # Include formula and exact mass in the output
            peak_lines.append(f"{mz:.6f} {intensity:.0f} \"{formula}\" {exact_mass:.6f}")
        else:
            peak_lines.append(f"{mz:.6f} {intensity:.0f}")
    
    msp_lines.append("\n".join(peak_lines))
    
    return "\n".join(msp_lines)

def write_output(msp_data: str, output_path: str) -> None:
    """
    Write the MSP data to a file.
    
    Args:
        msp_data: Formatted MSP string
        output_path: Path to write the MSP file
    """
    with open(output_path, 'w') as f:
        f.write(msp_data)
    logger.info(f"MSP file written to {output_path}")

def main():
    """Main function to run the MSP scraper."""
    parser = argparse.ArgumentParser(description='Scrape spectra from mzML files and convert to MSP format.')
    parser.add_argument('input_file', help='Path to the input TSV/Excel file.')
    parser.add_argument('--output_dir', default='output', help='Directory to store output MSP files.')
    parser.add_argument('--log_file', default='msp_scraper.log', help='Path to store log file.')
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Set up file handler for logging
    file_handler = logging.FileHandler(args.log_file)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(file_handler)
    
    logger.info(f"Processing input file: {args.input_file}")
    
    # Parse input file
    try:
        df = parse_input_file(args.input_file)
    except Exception as e:
        logger.error(f"Error parsing input file: {str(e)}")
        return 1
    
    # Process each row
    for idx, row in df.iterrows():
        try:
            # Check required fields
            if pd.isna(row.get('Molecule_name')):
                logger.error(f"Row {idx+1}: Missing required field 'Molecule_name'.")
                continue
            
            if pd.isna(row.get('mzML_filepath')):
                logger.error(f"Row {idx+1}: Missing required field 'mzML_filepath'.")
                continue
            
            # Extract spectrum
            spectrum_index = row.get('Spectrum_index') if pd.notna(row.get('Spectrum_index')) else None
            
            # Check for RT field - could be RT_seconds or RT
            retention_time = None
            if pd.notna(row.get('RT_seconds')):
                retention_time = float(row['RT_seconds'])
            elif pd.notna(row.get('RT')):
                # Handle RT in "min" format (e.g., "1.453 min")
                rt_str = str(row['RT'])
                if "min" in rt_str:
                    rt_str = rt_str.replace("min", "").strip()
                try:
                    # Convert minutes to seconds
                    retention_time = float(rt_str) * 60
                except ValueError:
                    logger.warning(f"Row {idx+1}: Could not parse RT value: {row['RT']}")
            
            ms_level = int(row['MS_level']) if pd.notna(row.get('MS_level')) else None
            
            if spectrum_index is None and retention_time is None:
                logger.error(f"Row {idx+1}: Either 'Spectrum_index' or retention time must be provided.")
                continue
            
            # Use index if available, otherwise use RT
            try:
                if spectrum_index is not None:
                    logger.info(f"Row {idx+1}: Extracting spectrum at index {spectrum_index}")
                    mz_array, intensity_array, spec_metadata = extract_spectrum(
                        row['mzML_filepath'], spectrum_index=int(spectrum_index), ms_level=ms_level
                    )
                else:
                    logger.info(f"Row {idx+1}: Finding spectrum closest to RT {retention_time} seconds")
                    mz_array, intensity_array, spec_metadata = extract_spectrum(
                        row['mzML_filepath'], retention_time=retention_time, ms_level=ms_level
                    )
            except Exception as e:
                logger.error(f"Row {idx+1}: Error extracting spectrum: {str(e)}")
                continue
            
            # Get intensity cutoff
            raw_intensity_cutoff = float(row['Raw_Intensity_Cutoff']) if pd.notna(row.get('Raw_Intensity_Cutoff')) else 0
            
            # Format MSP
            msp_data = format_msp(
                mz_array, intensity_array, spec_metadata, row, raw_intensity_cutoff
            )
            
            # Generate output filename
            molecule_name = str(row['Molecule_name']).replace(' ', '_')
            output_path = os.path.join(args.output_dir, f"{molecule_name}.msp")
            
            # Write output
            write_output(msp_data, output_path)
            
        except Exception as e:
            logger.error(f"Row {idx+1}: Error processing row - {str(e)}")
            logger.error(traceback.format_exc())
    
    logger.info("Processing complete.")
    return 0

if __name__ == "__main__":
    sys.exit(main()) 