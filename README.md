# MSPeeps

A Python package for extracting mass spectrometry spectra from mzML files and converting them to MSP format.

## Overview

This tool allows you to:

- Extract spectra from mzML files using either spectrum index or retention time
- Apply intensity cutoffs to filter peaks
- Convert SMILES to InChI and InChIKey when available
- Format the extracted data into MSP files according to standard conventions
- Process multiple spectra in batch mode via a tabular input file (TSV or Excel)
- Match peaks to molecular formulas using the fragfit package

## Installation

### From PyPI (Recommended)

Install the latest release from PyPI:

```bash
pip install mspeeps
```

### From conda

```bash
conda install -c gkreder mspeeps
```

### From Source

Clone the repository and install in development mode:

```bash
git clone https://github.com/gkreder/mspeeps.git
cd mspeeps
pip install -e .
```

### Using pixi

If you prefer using pixi for dependency management:

```bash
git clone https://github.com/gkreder/mspeeps.git
cd mspeeps
pixi install
```

## Usage

### Command-line Interface

Process an input file using the default settings:

```bash
mspeeps input_file.tsv
```

Specify custom output directory and log file:

```bash
mspeeps input_file.tsv --output_dir my_output --log_file custom_log.log
```

Enable verbose logging:

```bash
mspeeps input_file.tsv --verbose
```

### Python API

```python
import mspeeps
import pandas as pd

# Parse input file
df = mspeeps.parse_input_file("input_file.tsv")

# Process each row
for _, row in df.iterrows():
    output_path = mspeeps.process_file(row, output_dir="output")
    if output_path:
        print(f"Successfully processed: {output_path}")

# Or extract a spectrum directly
mz_array, intensity_array, metadata = mspeeps.extract_spectrum(
    mzml_path="file.mzML",
    spectrum_index=123,
    ms_level=2
)

# Format MSP data
msp_data = mspeeps.format_msp(
    mz_array,
    intensity_array,
    metadata,
    row_data={"Molecule_name": "Example"}
)
```

## Input Format

The input should be a TSV or Excel file with the following columns:

| Column Name | Description | Required? |
|-------------|-------------|-----------|
| Molecule_name | Name of the molecule | Yes |
| SMILES | SMILES notation | No |
| Molecular_formula | Chemical formula | No |
| Raw_Intensity_Cutoff | Cutoff for peak intensity | No (default: 0) |
| Formula_Matching_Tolerance | Tolerance for formula matching (in Da) | No |
| m/z | Precursor m/z | No |
| RT_seconds | Retention time in seconds | No* |
| RT | Retention time in minutes | No* |
| MS_level | MS level | No (default: 2) |
| Collision_energy | Collision energy used | No |
| mzML_filepath | Path to the mzML file | Yes |
| Spectrum_index | Index of the spectrum in the mzML file | No* |

\* Either Spectrum_index or retention time (RT_seconds/RT) must be provided.

Notes:
- If both spectrum index and RT are provided, the index is used.
- RT values in "min" format (e.g., "1.453 min") are automatically converted to seconds.

## Output Format

The output is an MSP file for each spectrum with the following format:

```
NAME: [Molecule_name]
[Additional metadata from input file]
INCHI: [Calculated from SMILES if provided]
INCHIKEY: [Calculated from SMILES if provided]
RETENTIONTIME: [Retention time in seconds]
PRECURSORMZ: [Precursor m/z]
MSLEVEL: [MS level]
NUM PEAKS: [Number of peaks]
[m/z] [intensity]
[m/z] [intensity]
...
```

## Formula Matching

The tool supports matching fragments in the spectrum to the closest possible molecular formula, within a specified tolerance, given the parent formula. This enables:

- **Fragment Formula Assignment**: Each m/z peak is annotated with its most likely molecular formula
- **Exact Mass Calculation**: The exact mass of each assigned formula is calculated
- **Enhanced Output Format**: Peak lines include formula, exact mass, and m/z difference (actual - theoretical): `[m/z] [intensity] "[formula]" [exact_mass] [m/z_difference]`

When formula matching is enabled, the output MSP file will look like:

```
NAME: Piperidine
SMILES: N1CCCCC1
MOLECULAR FORMULA: C5H11N
RAW INTENSITY CUTOFF: 100.0
FORMULA MATCHING TOLERANCE: 0.002
M/Z: 86.09643
RT SECONDS: 305.9
MS LEVEL: 2
COLLISION ENERGY: 0 V
MZML FILEPATH: /path/to/file.mzML
SPECTRUM INDEX: 576
INCHI: InChI=1S/C5H11N/c1-2-4-6-5-3-1/h6H,1-5H2
INCHIKEY: WEVYAHXRMPXWCK-UHFFFAOYSA-N
RETENTIONTIME: 5.04
PRECURSORMZ: 86.096430
MSLEVEL: 2
NUM PEAKS: 5
30.033819 2461 "CH4N" 30.033826 -0.000007
55.054611 1497 "C4H7" 55.054227 0.000384
57.070259 356 "C4H9" 57.069877 0.000382
68.049652 568 "C5H6" 68.049476 0.000176
84.080811 2834 "C5H10N" 84.080776 0.000035
```

## Development

### Running Tests

```bash
pytest
```

With coverage:

```bash
pytest --cov=mspeeps
```

## License

MIT

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

