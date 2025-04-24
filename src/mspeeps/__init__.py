"""
MSPeeps - Mass Spectrometry Peeps

A Python package for extracting mass spectrometry spectra from mzML files 
and converting them to MSP format.
"""

__version__ = "0.1.3"

# Import public API functions
from .core import extract_spectrum, format_msp, process_file
from .utils import convert_smiles_to_inchi, parse_input_file
from .formula_matching import match_formula

# Define what's available when using `from mspeeps import *`
__all__ = [
    "extract_spectrum",
    "format_msp", 
    "process_file",
    "convert_smiles_to_inchi",
    "parse_input_file",
    "match_formula",
] 