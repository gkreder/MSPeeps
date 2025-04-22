#!/usr/bin/env python3
"""
Test script for fragfit formula matching.

This script demonstrates how the fragfit package is used for formula matching in the MSP scraper.
"""

import numpy as np
import fragfit
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_formula_matching():
    """Test formula matching with a simple example."""
    # Example parent formula (Piperidine)
    parent_formula = "C5H11N"
    
    # Example m/z values (from the README example)
    mz_values = np.array([
        30.033819,
        55.054611,
        57.070259,
        68.049652,
        84.080811
    ])
    
    # Tolerance in Da
    tolerance = 0.002
    
    logger.info(f"Testing formula matching with parent formula: {parent_formula}")
    logger.info(f"Using tolerance: {tolerance} Da")
    
    results = []
    for mz in mz_values:
        # Use fragfit.find_best_form to find best matching formula
        formula, exact_mass, error = fragfit.find_best_form(
            mass=mz,
            parent_form=parent_formula,
            tolerance_da=tolerance,
            charge=1
        )
        
        results.append((mz, formula, exact_mass, error))
        logger.info(f"m/z: {mz:.6f} -> Formula: {formula}, Exact mass: {exact_mass:.6f}, Error: {error:.6f}")
    
    # Print results in MSP format
    print("\nMSP Format Example:")
    for mz, formula, exact_mass, _ in results:
        print(f"{mz:.6f} 1000 \"{formula}\" {exact_mass:.6f}")

if __name__ == "__main__":
    test_formula_matching() 