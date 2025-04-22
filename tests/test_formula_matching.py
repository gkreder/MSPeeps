"""
Tests for the formula matching functionality.
"""

import unittest
import numpy as np
from mspeeps.formula_matching import match_formula

class TestFormulaMatching(unittest.TestCase):
    """Test cases for formula matching functionality."""
    
    def test_match_formula(self):
        """Test the match_formula function with a simple example."""
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
        
        results = match_formula(mz_values, parent_formula, tolerance)
        
        # Check that we have the right number of results
        self.assertEqual(len(results), len(mz_values))
        
        # Check that each result has a formula and exact mass
        for i, (formula, exact_mass, error) in enumerate(results):
            self.assertIsNotNone(formula, f"Formula is None for m/z {mz_values[i]}")
            self.assertIsNotNone(exact_mass, f"Exact mass is None for m/z {mz_values[i]}")
            self.assertIsNotNone(error, f"Error is None for m/z {mz_values[i]}")
    
    def test_match_formula_invalid_input(self):
        """Test match_formula with invalid input."""
        # Test with invalid parent formula
        mz_values = np.array([100.0])
        results = match_formula(mz_values, "InvalidFormula", 0.002)
        
        # Should return a list with (None, None, None)
        self.assertEqual(len(results), 1)
        self.assertIsNone(results[0][0])
        self.assertIsNone(results[0][1])
        self.assertIsNone(results[0][2])

if __name__ == "__main__":
    unittest.main() 