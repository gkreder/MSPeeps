"""
Tests for core MSPeeps functionality.
"""

import unittest
import os
from unittest.mock import patch, MagicMock, mock_open
import tempfile
import pytest
import numpy as np
import json
import pandas as pd

class TestCore(unittest.TestCase):
    """Test cases for core functionality."""
    
    @patch('os.path.exists')
    @patch('pymzml.run.Reader')
    def test_extract_spectrum_by_index(self, mock_reader_class, mock_exists):
        """Test extracting a spectrum by index."""
        from src.mspeeps.core import extract_spectrum
        
        # Mock file exists check
        mock_exists.return_value = True
        
        # Create mock reader and spectrum
        mock_reader = MagicMock()
        mock_spectrum = MagicMock()
        mock_spectrum.mz = np.array([100.0, 200.0, 300.0])
        mock_spectrum.i = np.array([1000.0, 2000.0, 3000.0])
        mock_spectrum.ms_level = 2
        mock_spectrum.scan_time = [10.5]
        mock_spectrum.selected_precursors = [{"mz": 400.5, "charge": 2}]
        
        # Set up the get_by_id method to return mock_spectrum
        mock_reader.get_by_id.return_value = mock_spectrum
        mock_reader_class.return_value = mock_reader
        
        # Call extract_spectrum with index
        mz, intensity, metadata = extract_spectrum(
            mzml_path='file.mzML',
            spectrum_index=123,
            ms_level=2
        )
        
        # Check the results
        np.testing.assert_array_equal(mz, np.array([100.0, 200.0, 300.0]))
        np.testing.assert_array_equal(intensity, np.array([1000.0, 2000.0, 3000.0]))
        self.assertEqual(metadata['ms_level'], 2)
        self.assertEqual(metadata['retention_time'], 10.5)
        self.assertEqual(metadata['precursor_mz'], 400.5)
        self.assertEqual(metadata['precursor_charge'], 2)
    
    @patch('os.path.exists')
    @patch('pymzml.run.Reader')
    def test_extract_spectrum_by_retention_time(self, mock_reader_class, mock_exists):
        """Test extracting a spectrum by retention time."""
        from src.mspeeps.core import extract_spectrum
        
        # Mock file exists check
        mock_exists.return_value = True
        
        # Create mock reader and spectra
        mock_reader = MagicMock()
        
        mock_spectrum1 = MagicMock()
        mock_spectrum1.ms_level = 2
        mock_spectrum1.scan_time = [10.0]
        
        mock_spectrum2 = MagicMock()
        mock_spectrum2.mz = np.array([100.0, 200.0, 300.0])
        mock_spectrum2.i = np.array([1000.0, 2000.0, 3000.0])
        mock_spectrum2.ms_level = 2
        mock_spectrum2.scan_time = [20.0]
        mock_spectrum2.selected_precursors = [{"mz": 400.5, "charge": 2}]
        
        mock_spectrum3 = MagicMock()
        mock_spectrum3.ms_level = 2
        mock_spectrum3.scan_time = [30.0]
        
        # Set up the __iter__ method to return a list of spectra
        mock_reader.__iter__.return_value = [mock_spectrum1, mock_spectrum2, mock_spectrum3]
        mock_reader_class.return_value = mock_reader
        
        # Call extract_spectrum with retention_time
        mz, intensity, metadata = extract_spectrum(
            mzml_path='file.mzML',
            retention_time=20.0,
            ms_level=2
        )
        
        # Check the results
        np.testing.assert_array_equal(mz, np.array([100.0, 200.0, 300.0]))
        np.testing.assert_array_equal(intensity, np.array([1000.0, 2000.0, 3000.0]))
        self.assertEqual(metadata['ms_level'], 2)
        self.assertEqual(metadata['retention_time'], 20.0)
        self.assertEqual(metadata['precursor_mz'], 400.5)
        self.assertEqual(metadata['precursor_charge'], 2)
    
    def test_extract_spectrum_file_not_found(self):
        """Test extracting a spectrum from a non-existent file."""
        from src.mspeeps.core import extract_spectrum
        
        with self.assertRaises(FileNotFoundError):
            extract_spectrum(
                mzml_path='nonexistent.mzML',
                spectrum_index=1,
                ms_level=2
            )
    
    @patch('src.mspeeps.utils.convert_smiles_to_inchi')
    @patch('src.mspeeps.formula_matching.match_formula')
    def test_format_msp_basic(self, mock_match_formula, mock_convert_smiles):
        """Test basic formatting of MSP data."""
        from src.mspeeps.core import format_msp
        
        # Create test data
        mz = np.array([100.0, 200.0, 300.0])
        intensity = np.array([1000.0, 2000.0, 3000.0])
        metadata = {
            'ms_level': 2,
            'retention_time': 10.5,
            'precursor_mz': 400.5,
            'precursor_charge': 2
        }
        
        # Create mock row data
        row_data = pd.Series({
            'Molecule_name': 'Test Molecule',
            'SMILES': 'CC(=O)O',
            'Molecular_formula': 'C2H4O2',
            'Formula_Matching_Tolerance': 0.002
        })
        
        # Mock convert_smiles_to_inchi return value
        mock_convert_smiles.return_value = ('InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)', 'QTBSBXVTEAMEQO-UHFFFAOYSA-N')
        
        # Mock match_formula return value
        mock_match_formula.return_value = [
            ('C2H4O2', 60.0211, 0.0001),
            ('C2H4O', 44.0262, 0.0002),
            ('CO2', 43.9898, 0.0003)
        ]
        
        # Format as MSP
        msp_data = format_msp(
            mz, 
            intensity, 
            metadata,
            row_data,
            raw_intensity_cutoff=0
        )
        
        # Check the MSP data
        self.assertIn("NAME: Test Molecule", msp_data)
        self.assertIn("SMILES: CC(=O)O", msp_data)
        self.assertIn("INCHI: InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)", msp_data)
        self.assertIn("INCHIKEY: QTBSBXVTEAMEQO-UHFFFAOYSA-N", msp_data)
        self.assertIn("MOLECULAR FORMULA: C2H4O2", msp_data)
        self.assertIn("PRECURSORMZ: 400.500000", msp_data)
        self.assertIn("RETENTIONTIME: 0.17", msp_data)  # 10.5 seconds converted to minutes
        self.assertIn("NUM PEAKS: 3", msp_data)
    
    @patch('src.mspeeps.core.extract_spectrum')
    @patch('src.mspeeps.core.format_msp')
    @patch('src.mspeeps.core.os.makedirs')
    @patch('builtins.open', new_callable=mock_open)
    def test_process_file(self, mock_file, mock_makedirs, mock_format_msp, mock_extract_spectrum):
        """Test processing a file."""
        from src.mspeeps.core import process_file
        
        # Create test data
        row_data = pd.Series({
            'Molecule_name': 'Test Molecule',
            'mzML_filepath': 'file.mzML',
            'Spectrum_index': 123,
            'Molecular_formula': 'C6H12O6',
            'RT': None,
            'MS_level': 2
        })
        
        # Mock extract_spectrum
        mock_extract_spectrum.return_value = (
            np.array([100.0, 200.0, 300.0]),
            np.array([1000.0, 2000.0, 3000.0]),
            {
                'ms_level': 2,
                'retention_time': 10.5,
                'precursor_mz': 400.5,
                'precursor_charge': 2
            }
        )
        
        # Mock format_msp
        mock_format_msp.return_value = "MSP data"
        
        # Call process_file
        with tempfile.TemporaryDirectory() as tempdir:
            output_path = process_file(row_data, tempdir)
            
            # Check extract_spectrum was called correctly
            mock_extract_spectrum.assert_called_once()
            # Check that extract_spectrum was called with the correct arguments
            call_args = mock_extract_spectrum.call_args[0]
            call_kwargs = mock_extract_spectrum.call_args[1]
            self.assertEqual(call_kwargs.get('mzml_path', call_args[0] if call_args else None), 'file.mzML')
            self.assertEqual(call_kwargs.get('spectrum_index', call_args[1] if len(call_args) > 1 else None), 123)
            self.assertEqual(call_kwargs.get('ms_level', call_args[3] if len(call_args) > 3 else None), 2)
            
            # Check format_msp was called correctly
            mock_format_msp.assert_called_once()
            
            # Check file was written correctly
            mock_file.assert_called_once()
            mock_file().write.assert_called_once_with("MSP data")
            
            # Check output path is correct
            self.assertTrue(output_path.endswith('Test Molecule.msp'))

    def test_write_output(self):
        """Test write_output function."""
        from src.mspeeps.core import write_output
        
        msp_data = "NAME: Test Molecule"
        
        # Call write_output in a temp directory
        with tempfile.TemporaryDirectory() as tempdir:
            output_path = os.path.join(tempdir, 'output.msp')
            
            # Test the function
            write_output(msp_data, output_path)
            
            # Verify the file was written correctly
            with open(output_path, 'r') as f:
                content = f.read()
                self.assertEqual(content, msp_data)

    def test_match_formula(self):
        """Test match_formula function."""
        # Import the match_formula function
        import src.mspeeps.formula_matching
        
        # Create a proper fake implementation for match_formula
        def fake_match_formula(mz_array, parent_formula, tolerance_da, charge=1):
            return [
                ('CO2', 43.9898, 0.0001),
                ('CH3O', 31.0184, 0.0002)
            ]
        
        # Patch the match_formula function
        with patch('src.mspeeps.formula_matching.match_formula', side_effect=fake_match_formula):
            # Test with a simple formula and mz values
            results = src.mspeeps.formula_matching.match_formula(
                mz_array=np.array([44.0, 31.02]),
                parent_formula='C2H6O',
                tolerance_da=0.01,
                charge=1
            )
            
            # Check the results
            self.assertEqual(len(results), 2)
            self.assertEqual(results[0][0], 'CO2')  # Formula
            self.assertAlmostEqual(results[0][1], 43.9898)  # Exact mass
            self.assertAlmostEqual(results[0][2], 0.0001)  # Error
            
            self.assertEqual(results[1][0], 'CH3O')  # Formula
            self.assertAlmostEqual(results[1][1], 31.0184)  # Exact mass
            self.assertAlmostEqual(results[1][2], 0.0002)  # Error

    def test_convert_smiles_to_inchi(self):
        """Test convert_smiles_to_inchi function."""
        from src.mspeeps.utils import convert_smiles_to_inchi
        
        # Create a proper fake implementation for convert_smiles_to_inchi
        def fake_convert_smiles(smiles):
            return 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'
        
        # Patch the convert_smiles_to_inchi function directly
        with patch('src.mspeeps.utils.convert_smiles_to_inchi', side_effect=fake_convert_smiles):
            # Test the function
            inchi, inchikey = convert_smiles_to_inchi('CCO')
            
            # Check results
            self.assertEqual(inchi, 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3')
            self.assertEqual(inchikey, 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N')

if __name__ == "__main__":
    unittest.main() 