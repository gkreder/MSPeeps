"""
Tests for CLI interface.
"""

import unittest
from unittest.mock import patch, MagicMock, call
import tempfile
import os
import sys
import numpy as np
import pandas as pd
from io import StringIO

from src.mspeeps.cli import (
    main,
    batch_command,
    extract_command,
    convert_command,
    match_formula_command,
    convert_smiles_command,
    info_command
)

class TestCLI(unittest.TestCase):
    """Test cases for CLI interface."""
    
    def test_handle_version(self):
        """Test version command output."""
        # Redirect stdout to capture output
        captured_output = StringIO()
        sys.stdout = captured_output
        
        # Patch the version module
        with patch('src.mspeeps.cli.setup_logging'), \
             patch('src.mspeeps.__version__', 'test_version'):
            # Call the version handler
            main(['--version'])
        
        # Reset stdout
        sys.stdout = sys.__stdout__
        
        # Check that the output contains version information
        output = captured_output.getvalue()
        self.assertIn('MSPeeps', output)
        self.assertIn('version', output)
    
    @patch('src.mspeeps.cli.parse_input_file')
    @patch('src.mspeeps.cli.process_file')
    @patch('src.mspeeps.cli.os.makedirs')
    @patch('src.mspeeps.cli.os.path.exists')
    def test_batch_command(self, mock_exists, mock_makedirs, mock_process_file, mock_parse_input_file):
        """Test batch command."""
        # Create mock DataFrame
        mock_df = pd.DataFrame({
            'Molecule_name': ['Compound1', 'Compound2'],
            'mzML_filepath': ['file1.mzML', 'file2.mzML'],
            'Spectrum_index': [1, 2],
            'Precursor_mz': [100.0, 200.0],
            'Precursor_charge': [1, 2]
        })
        mock_parse_input_file.return_value = mock_df
        
        # Configure mock return values
        mock_exists.return_value = False  # To trigger makedirs
        mock_process_file.return_value = 'output_dir/test.msp'
            
        # Create a mock logger
        mock_logger = MagicMock()
        
        # Create mock args
        args = MagicMock()
        args.input_file = 'test_input.csv'
        args.output_dir = 'output_dir'
        args.verbose = False
        
        # Call the batch command
        batch_command(args, mock_logger)
        
        # Check parse_input_file was called
        mock_parse_input_file.assert_called_once_with('test_input.csv')
        
        # Check the output directory was created
        mock_makedirs.assert_called_once_with('output_dir')
        
        # Check process_file was called for each row
        self.assertEqual(mock_process_file.call_count, 2)
    
    @patch('src.mspeeps.cli.extract_spectrum')
    @patch('src.mspeeps.cli.os.path.exists')
    def test_extract_command(self, mock_exists, mock_extract_spectrum):
        """Test extract command."""
        # Setup mocks
        mock_exists.return_value = True
        mock_extract_spectrum.return_value = (
            np.array([100.0, 200.0, 300.0]),           # mz
            np.array([1000.0, 2000.0, 3000.0]),        # intensity
            {                                          # metadata
                'ms_level': 2,
                'retention_time': 100.0,
                'precursor_mz': 400.0,
                'precursor_charge': 1
            }
        )
        
        # Create a mock logger
        mock_logger = MagicMock()
        
        # Create mock args
        args = MagicMock()
        args.mzml_file = 'test.mzML'
        args.spectrum_index = 0
        args.retention_time = None
        args.ms_level = 2
        args.output = 'output.msp'
        args.format = 'json'
        
        # Set up open mock
        with patch('builtins.open', unittest.mock.mock_open()) as mock_open:
            # Call the extract command
            extract_command(args, mock_logger)
            
            # Check extract_spectrum was called with the correct parameters
            mock_extract_spectrum.assert_called_once_with(
                mzml_path='test.mzML', 
                spectrum_index=0, 
                retention_time=None, 
                ms_level=2
            )
    
    @patch('src.mspeeps.cli.format_msp')
    @patch('src.mspeeps.cli.np.array')
    @patch('src.mspeeps.cli.os.path.exists')
    def test_convert_command(self, mock_exists, mock_array, mock_format_msp):
        """Test convert command."""
        # Setup mocks
        mock_exists.return_value = True
        mock_array.side_effect = [
            np.array([100.0, 200.0, 300.0]),
            np.array([1000.0, 2000.0, 3000.0])
        ]
        mock_format_msp.return_value = "NAME: Test Compound\n"
        
        # Create a mock logger
        mock_logger = MagicMock()
        
        # Create mock args
        args = MagicMock()
        args.json_file = None
        args.mz_file = 'test_mz.txt'
        args.intensity_file = 'test_intensity.txt'
        args.metadata_file = None
        args.output = 'output.msp'
        args.intensity_cutoff = 0
        
        # Set up open mock with mock data
        mock_mz_data = "100.0\n200.0\n300.0"
        mock_intensity_data = "1000.0\n2000.0\n3000.0"
        
        # Mock open to handle sequential file reads
        open_mock = unittest.mock.mock_open()
        handles = [
            unittest.mock.mock_open(read_data=mock_mz_data).return_value,
            unittest.mock.mock_open(read_data=mock_intensity_data).return_value,
            unittest.mock.mock_open().return_value  # For output
        ]
        open_mock.side_effect = handles
            
        # Call the convert command with mocked open
        with patch('builtins.open', open_mock):
            convert_command(args, mock_logger)
            
            # Check format_msp was called
            mock_format_msp.assert_called_once()
    
    @patch('src.mspeeps.cli.match_formula')
    def test_match_formula_command(self, mock_match_formula):
        """Test match_formula command."""
        # Setup mock for match_formula
        mock_match_formula.return_value = [
            ('C2H4', 28.0313, 0.0002)
        ]
        
        # Create a mock logger
        mock_logger = MagicMock()
        
        # Create mock args
        args = MagicMock()
        args.mz_values = '28.0315'  # As a string
        args.parent_formula = 'C2H6O'
        args.tolerance = 0.001
        args.charge = 1
        args.output = None
        args.format = 'tsv'
        
        # Mock np.array to handle string-to-array conversion
        with patch('src.mspeeps.cli.np.array', return_value=np.array([28.0315])):
            # Redirect stdout to capture output
            captured_output = StringIO()
            sys.stdout = captured_output
            
            # Call the match formula command
            match_formula_command(args, mock_logger)
            
            # Reset stdout
            sys.stdout = sys.__stdout__
            
            # Check match_formula was called
            mock_match_formula.assert_called_once()
    
    @patch('src.mspeeps.cli.convert_smiles_to_inchi')
    def test_convert_smiles_command(self, mock_convert_to_inchi):
        """Test convert_smiles command."""
        # Setup mocks
        mock_convert_to_inchi.return_value = ('InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N')
        
        # Create a mock logger
        mock_logger = MagicMock()
        
        # Create mock args
        args = MagicMock()
        args.smiles = 'CCO'
        args.output = None
        
        # Redirect stdout to capture output
        captured_output = StringIO()
        sys.stdout = captured_output
        
        # Call the convert_smiles command
        convert_smiles_command(args, mock_logger)
        
        # Reset stdout
        sys.stdout = sys.__stdout__
        
        # Check the mocks were called correctly
        mock_convert_to_inchi.assert_called_once_with('CCO')
        
        # Check output contains the correct information
        output = captured_output.getvalue()
        self.assertIn('CCO', output)
        self.assertIn('InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', output)
        self.assertIn('LFQSCWFLJHTTHZ-UHFFFAOYSA-N', output)
    
    def test_info_command(self):
        """Test info command."""
        # Create mock spectra
        mock_spectrum1 = MagicMock()
        mock_spectrum1.ms_level = 1
        mock_spectrum1.scan_time = (100.0,)
        
        mock_spectrum2 = MagicMock()
        mock_spectrum2.ms_level = 2
        mock_spectrum2.scan_time = (200.0,)
        
        mock_spectrum3 = MagicMock()
        mock_spectrum3.ms_level = 2
        mock_spectrum3.scan_time = (300.0,)
        
        # Setup mock for pymzml.run.Reader
        mock_reader_instance = MagicMock()
        mock_reader_instance.__iter__.return_value = [mock_spectrum1, mock_spectrum2, mock_spectrum3]
        
        # Create a mock logger
        mock_logger = MagicMock()
        
        # Create mock args
        args = MagicMock()
        args.mzml_file = 'test.mzML'
        args.output = None
        args.format = 'text'
        
        # Mock os.path.exists and os.path.getsize
        with patch('os.path.exists', return_value=True), \
             patch('os.path.getsize', return_value=1024*1024), \
             patch('pymzml.run.Reader', return_value=mock_reader_instance):
            
            # Redirect stdout to capture output
            captured_output = StringIO()
            sys.stdout = captured_output
            
            # Call the info command
            info_command(args, mock_logger)
            
            # Reset stdout
            sys.stdout = sys.__stdout__
            
            # Check output contains the correct information
            output = captured_output.getvalue()
            self.assertIn('test.mzML', output)
            self.assertIn('Spectra: 3', output)
            self.assertIn('MS1 Spectra: 1', output)
            self.assertIn('MS2 Spectra: 2', output)

if __name__ == "__main__":
    unittest.main() 