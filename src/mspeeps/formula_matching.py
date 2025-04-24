"""
Formula matching functionality for MS spectra.

This module handles matching of m/z values to molecular formulas.
"""

import logging
import numpy as np
from typing import List, Tuple, Optional, Union

# Set up logging
logger = logging.getLogger(__name__)

def match_formula(mz_array: np.ndarray, parent_formula: str, 
                 tolerance_da: float, charge: int = 1) -> List[Tuple[Optional[str], Optional[float], Optional[float]]]:
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
    try:
        import fragfit
    except ImportError:
        logger.error("fragfit package not found. Install with `pip install fragfit`")
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

def test_formula_matching(parent_formula: str, mz_values: np.ndarray, tolerance: float = 0.002) -> None:
    """
    Test formula matching with provided values.
    
    Args:
        parent_formula: Molecular formula of the parent molecule
        mz_values: Array of m/z values to match
        tolerance: Mass tolerance in Da (default: 0.002)
    """
    logger.info(f"Testing formula matching with parent formula: {parent_formula}")
    logger.info(f"Using tolerance: {tolerance} Da")
    
    results = match_formula(mz_values, parent_formula, tolerance)
    
    # Print results 
    for mz, (formula, exact_mass, error) in zip(mz_values, results):
        if formula and exact_mass:
            logger.info(f"m/z: {mz:.6f} -> Formula: {formula}, Exact mass: {exact_mass:.6f}, Error: {error:.6f}")
        else:
            logger.info(f"m/z: {mz:.6f} -> No match found within tolerance")
    
    # Print results in MSP format
    print("\nMSP Format Example:")
    for mz, (formula, exact_mass, error) in zip(mz_values, results):
        if formula and exact_mass:
            mz_difference = mz - exact_mass
            print(f"{mz:.6f} 1000 \"{formula}\" {exact_mass:.6f} {mz_difference:.6f}")
        else:
            print(f"{mz:.6f} 1000") 