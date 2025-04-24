"""
Core functionality for extracting and processing mass spectrometry spectra.
"""

import os
import numpy as np
import logging
import pandas as pd
from typing import Dict, List, Tuple, Optional, Union, Any

from .utils import convert_smiles_to_inchi, parse_input_file
from .formula_matching import match_formula

# Set up logging
logger = logging.getLogger(__name__)

def extract_spectrum(mzml_path: str, spectrum_index: Optional[int] = None, 
                    retention_time: Optional[float] = None, 
                    ms_level: int = 2) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
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

def format_msp(mz_array: np.ndarray, intensity_array: np.ndarray, 
              metadata: Dict[str, Any], row_data: pd.Series, 
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
            logger.error(f"Error in formula matching: {str(e)}")
            formula_matches = None
    
    # Create MSP string
    msp_lines = []
    
    # Add metadata
    msp_lines.append(f"NAME: {row_data.get('Molecule_name', 'Unknown')}")
    
    # Add SMILES if available
    if pd.notna(row_data.get('SMILES')):
        msp_lines.append(f"SMILES: {row_data['SMILES']}")
        
        # Add InChI and InChIKey if SMILES is available
        inchi, inchikey = convert_smiles_to_inchi(row_data['SMILES'])
        if inchi:
            msp_lines.append(f"INCHI: {inchi}")
            msp_lines.append(f"INCHIKEY: {inchikey}")
    
    # Add molecular formula if available
    if pd.notna(row_data.get('Molecular_formula')):
        msp_lines.append(f"MOLECULAR FORMULA: {row_data['Molecular_formula']}")
    
    # Add other metadata from the input file
    if pd.notna(row_data.get('Raw_Intensity_Cutoff')):
        msp_lines.append(f"RAW INTENSITY CUTOFF: {row_data['Raw_Intensity_Cutoff']}")
    
    if pd.notna(row_data.get('Formula_Matching_Tolerance')):
        msp_lines.append(f"FORMULA MATCHING TOLERANCE: {row_data['Formula_Matching_Tolerance']}")
    
    if pd.notna(row_data.get('m/z')):
        msp_lines.append(f"M/Z: {row_data['m/z']}")
    
    # Convert RT from minutes to seconds if needed and add to MSP
    rt_seconds = None
    if pd.notna(row_data.get('RT_seconds')):
        rt_seconds = row_data['RT_seconds']
    elif pd.notna(row_data.get('RT')):
        # Check if the RT is a string with "min" suffix
        rt = row_data['RT']
        if isinstance(rt, str) and "min" in rt:
            try:
                rt = float(rt.strip().split()[0])
            except (ValueError, IndexError):
                pass
        rt_seconds = float(rt) * 60  # Convert minutes to seconds
    
    if rt_seconds is not None:
        msp_lines.append(f"RT SECONDS: {rt_seconds}")
    
    # Add metadata from the spectrum
    if metadata.get('ms_level') is not None:
        msp_lines.append(f"MS LEVEL: {metadata['ms_level']}")
    
    # Add collision energy if available
    if pd.notna(row_data.get('Collision_energy')):
        msp_lines.append(f"COLLISION ENERGY: {row_data['Collision_energy']}")
    
    # Add mzML filepath and spectrum index
    if pd.notna(row_data.get('mzML_filepath')):
        msp_lines.append(f"MZML FILEPATH: {row_data['mzML_filepath']}")
    
    if pd.notna(row_data.get('Spectrum_index')):
        msp_lines.append(f"SPECTRUM INDEX: {row_data['Spectrum_index']}")
    
    # Add retention time info
    rt = metadata.get('retention_time')
    if rt is not None:
        rt_min = rt / 60.0  # Convert seconds to minutes
        msp_lines.append(f"RETENTIONTIME: {rt_min:.2f}")
    
    # Add precursor m/z
    precursor_mz = metadata.get('precursor_mz')
    if precursor_mz is not None:
        msp_lines.append(f"PRECURSORMZ: {precursor_mz:.6f}")
    
    # Add MS level again (required by some software)
    msp_lines.append(f"MSLEVEL: {metadata.get('ms_level', 2)}")
    
    # Add peak count
    msp_lines.append(f"NUM PEAKS: {len(mz_array)}")
    
    # Add peaks
    if formula_matches and do_formula_matching:
        for i, (mz, intensity) in enumerate(zip(mz_array, intensity_array)):
            formula, exact_mass, _ = formula_matches[i]
            if formula and exact_mass:
                # Calculate and include the m/z difference (actual - theoretical)
                mz_difference = mz - exact_mass
                msp_lines.append(f"{mz:.6f} {intensity:.0f} \"{formula}\" {exact_mass:.6f} {mz_difference:.6f}")
            else:
                msp_lines.append(f"{mz:.6f} {intensity:.0f}")
    else:
        for mz, intensity in zip(mz_array, intensity_array):
            msp_lines.append(f"{mz:.6f} {intensity:.0f}")
    
    return "\n".join(msp_lines)

def write_output(msp_data: str, output_path: str) -> None:
    """
    Write MSP data to output file.
    
    Args:
        msp_data: Formatted MSP data
        output_path: Path to write output file
    """
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Write to file
    with open(output_path, 'w') as f:
        f.write(msp_data)
    
    logger.info(f"Wrote MSP data to {output_path}")

def process_file(row_data: pd.Series, output_dir: str = "output") -> Optional[str]:
    """
    Process a single row from the input file and create MSP output.
    
    Args:
        row_data: Row data from the input file
        output_dir: Directory to store output files

    Returns:
        Path to the output file if successful, None otherwise
    """
    try:
        # Get required fields
        molecule_name = row_data.get('Molecule_name')
        mzml_filepath = row_data.get('mzML_filepath')
        
        if pd.isna(molecule_name) or pd.isna(mzml_filepath):
            logger.error("Missing required fields: Molecule_name or mzML_filepath")
            return None
        
        # Get optional fields with defaults
        spectrum_index = row_data.get('Spectrum_index')
        if pd.notna(spectrum_index):
            spectrum_index = int(spectrum_index)
        else:
            spectrum_index = None
        
        # Get retention time
        rt_seconds = None
        if pd.notna(row_data.get('RT_seconds')):
            rt_seconds = float(row_data['RT_seconds'])
        elif pd.notna(row_data.get('RT')):
            # Check if the RT is a string with "min" suffix
            rt = row_data['RT']
            if isinstance(rt, str) and "min" in rt:
                try:
                    rt = float(rt.strip().split()[0])
                except (ValueError, IndexError):
                    pass
            rt_seconds = float(rt) * 60  # Convert minutes to seconds
        
        # Get MS level
        ms_level = row_data.get('MS_level', 2)
        if pd.isna(ms_level):
            ms_level = 2
        else:
            ms_level = int(ms_level)
        
        # Get intensity cutoff
        raw_intensity_cutoff = row_data.get('Raw_Intensity_Cutoff', 0)
        if pd.isna(raw_intensity_cutoff):
            raw_intensity_cutoff = 0
        else:
            raw_intensity_cutoff = float(raw_intensity_cutoff)
        
        # Extract spectrum
        logger.info(f"Extracting spectrum from {mzml_filepath}")
        if spectrum_index is not None:
            logger.info(f"Using spectrum index: {spectrum_index}")
            mz_array, intensity_array, metadata = extract_spectrum(
                mzml_filepath, 
                spectrum_index=spectrum_index, 
                ms_level=ms_level
            )
        elif rt_seconds is not None:
            logger.info(f"Using retention time: {rt_seconds} seconds")
            mz_array, intensity_array, metadata = extract_spectrum(
                mzml_filepath, 
                retention_time=rt_seconds, 
                ms_level=ms_level
            )
        else:
            logger.error("Either Spectrum_index or retention time (RT_seconds/RT) must be provided")
            return None
        
        # Format MSP
        logger.info("Formatting MSP data")
        msp_data = format_msp(
            mz_array, 
            intensity_array, 
            metadata, 
            row_data, 
            raw_intensity_cutoff
        )
        
        # Write output
        output_filename = f"{molecule_name}.msp"
        output_path = os.path.join(output_dir, output_filename)
        write_output(msp_data, output_path)
        
        return output_path
    
    except Exception as e:
        logger.error(f"Error processing row: {str(e)}")
        return None 