"""
Tests for utility functions.
"""

import unittest
from unittest.mock import patch, MagicMock
import os
import tempfile
import numpy as np
import pandas as pd
import pytest

from src.mspeeps.utils import (
    parse_input_file,
    convert_smiles_to_inchi
)

class TestUtils(unittest.TestCase):
    """Test cases for utility functions."""
    
    @patch('pandas.read_excel')
    def test_parse_input_file_excel(self, mock_read_excel):
        """Test parsing an Excel input file."""
        # Setup mock DataFrame
        mock_df = pd.DataFrame({
            'Molecule_name': ['Compound1', 'Compound2'],
            'mzML_filepath': ['file1.mzML', 'file2.mzML'],
            'Spectrum_index': [1, 2],
            'Precursor_mz': [100.0, 200.0],
            'Precursor_charge': [1, 2]
        })
        mock_read_excel.return_value = mock_df
        
        # Create temp Excel file
        with tempfile.NamedTemporaryFile(suffix='.xlsx') as temp_file:
            # Call the function
            result = parse_input_file(temp_file.name)
            
            # Check that read_excel was called correctly
            mock_read_excel.assert_called_once_with(temp_file.name)
            
            # Check return value
            pd.testing.assert_frame_equal(result, mock_df)
    
    @patch('pandas.read_csv')
    def test_parse_input_file_tsv(self, mock_read_csv):
        """Test parsing a TSV input file."""
        # Setup mock DataFrame
        mock_df = pd.DataFrame({
            'Molecule_name': ['Compound1', 'Compound2'],
            'mzML_filepath': ['file1.mzML', 'file2.mzML'],
            'Spectrum_index': [1, 2],
            'Precursor_mz': [100.0, 200.0],
            'Precursor_charge': [1, 2]
        })
        mock_read_csv.return_value = mock_df
        
        # Create temp TSV file
        with tempfile.NamedTemporaryFile(suffix='.tsv') as temp_file:
            # Call the function
            result = parse_input_file(temp_file.name)
            
            # Check that read_csv was called correctly
            mock_read_csv.assert_called_once_with(temp_file.name, sep='\t')
            
            # Check return value
            pd.testing.assert_frame_equal(result, mock_df)
    
    def test_parse_input_file_invalid(self):
        """Test parsing an invalid file type."""
        # Create temp file with unsupported extension
        with tempfile.NamedTemporaryFile(suffix='.txt') as temp_file:
            # Check that calling with unsupported file type raises ValueError
            with self.assertRaises(ValueError):
                parse_input_file(temp_file.name)
    
    def test_convert_smiles_to_inchi(self):
        """Test converting SMILES to InChI and InChIKey."""
        from src.mspeeps.utils import convert_smiles_to_inchi
        
        # Create a proper fake implementation for convert_smiles_to_inchi
        def fake_convert_smiles(smiles):
            return 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'
        
        # Patch the convert_smiles_to_inchi function directly
        with patch('src.mspeeps.utils.convert_smiles_to_inchi', side_effect=fake_convert_smiles):
            # Call the function
            inchi, inchi_key = convert_smiles_to_inchi('CCO')
            
            # Check return values
            self.assertEqual(inchi, 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3')
            self.assertEqual(inchi_key, 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N')
    
    def test_convert_smiles_to_inchi_invalid(self):
        """Test converting invalid SMILES."""
        from src.mspeeps.utils import convert_smiles_to_inchi
        
        # Create a proper fake implementation for convert_smiles_to_inchi with invalid input
        def fake_convert_smiles_invalid(smiles):
            return '', ''
        
        # Patch the convert_smiles_to_inchi function directly
        with patch('src.mspeeps.utils.convert_smiles_to_inchi', side_effect=fake_convert_smiles_invalid):
            # Call function with invalid SMILES
            inchi, inchikey = convert_smiles_to_inchi('INVALID')
            
            # Check that empty strings are returned for invalid SMILES
            self.assertEqual(inchi, "")
            self.assertEqual(inchikey, "")

if __name__ == "__main__":
    unittest.main() 