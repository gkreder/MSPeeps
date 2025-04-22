"""
Utility functions for MSPeeps.

This module contains helper functions for file parsing, SMILES conversion, etc.
"""

import logging
import pandas as pd
from typing import Tuple, Dict, Any

# Set up logging
logger = logging.getLogger(__name__)

def parse_input_file(file_path: str) -> pd.DataFrame:
    """
    Parse the input TSV/Excel file.
    
    Args:
        file_path: Path to the input file (TSV or Excel)
    
    Returns:
        DataFrame containing the parsed data
    """
    if file_path.endswith('.tsv'):
        df = pd.read_csv(file_path, sep='\t')
    elif file_path.endswith(('.xlsx', '.xls')):
        df = pd.read_excel(file_path)
    else:
        raise ValueError("Input file must be a TSV or Excel file.")
    
    # Check required columns
    required_cols = ['Molecule_name', 'mzML_filepath']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {', '.join(missing_cols)}")
    
    # Check that each row has either Spectrum_index or RT/RT_seconds
    has_spectrum_index = df['Spectrum_index'].notna()
    has_rt = df['RT'].notna() if 'RT' in df.columns else pd.Series(False, index=df.index)
    has_rt_seconds = df['RT_seconds'].notna() if 'RT_seconds' in df.columns else pd.Series(False, index=df.index)
    
    missing_index_and_rt = ~(has_spectrum_index | has_rt | has_rt_seconds)
    if missing_index_and_rt.any():
        missing_rows = df[missing_index_and_rt].index.tolist()
        logger.warning(f"Rows {missing_rows} are missing both Spectrum_index and RT/RT_seconds")
    
    return df

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
    except ImportError:
        logger.error("RDKit not found. Install with `pip install rdkit`")
        return "", ""
    
    try:
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

def setup_logging(log_file: str = None, log_level: int = logging.INFO) -> logging.Logger:
    """
    Set up logging configuration.
    
    Args:
        log_file: Path to log file (if None, logs to console only)
        log_level: Logging level
    
    Returns:
        Configured logger
    """
    logger = logging.getLogger("mspeeps")
    logger.setLevel(log_level)
    
    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # Create console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # Create file handler if log_file is specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger 