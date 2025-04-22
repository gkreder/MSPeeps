# MSP File Scraper

A Python package for extracting mass spectrometry spectra from mzML files and converting them to MSP format.

## Overview

This tool allows you to:

- Extract spectra from mzML files using either spectrum index or retention time
- Apply intensity cutoffs to filter peaks
- Convert SMILES to InChI and InChIKey when available
- Format the extracted data into MSP files according to standard conventions
- Process multiple spectra in batch mode via a tabular input file (TSV or Excel)

## Installation

This package uses [pixi](https://prefix.dev/) for dependency management. Make sure you have pixi installed before proceeding.

### Clone the repository:

```bash
git clone https://github.com/yourusername/MSP-file-scraping.git
cd MSP-file-scraping
```

### Install dependencies:

```bash
pixi install
```

Or install dependencies individually:

```bash
pixi add pandas
pixi add --pypi pymzml
pixi add --pypi rdkit
pixi add --pypi openpyxl
```

## Dependencies

- **pandas**: For reading and processing tabular data
- **pymzml**: For parsing mzML files
- **RDKit**: For SMILES to InChI/InChIKey conversion
- **openpyxl**: For Excel file support
- **numpy**: For numerical operations

## Usage

### Basic Usage

Process an input file using the default settings:

```bash
pixi run --manifest-path /path/to/pixi.toml python msp_scraper.py input_file.tsv
```

### With Custom Settings

Specify custom output directory and log file:

```bash
pixi run --manifest-path /path/to/pixi.toml python msp_scraper.py input_file.tsv --output_dir my_output --log_file custom_log.log
```

### Command-line Arguments

- `input_file`: Path to the input TSV/Excel file (required)
- `--output_dir`: Directory to store output MSP files (default: 'output')
- `--log_file`: Path to store log file (default: 'msp_scraper.log')

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

### Future Enhancement: Formula Matching

In future versions, the tool will support matching fragments in the spectrum to the closest possible molecular formula, within a specified tolerance, given the parent formula. This will enable:

- **Fragment Formula Assignment**: Each m/z peak will be annotated with its most likely molecular formula
- **Exact Mass Calculation**: The exact mass of each assigned formula will be calculated
- **Format Enhancement**: Peak lines will include formula and exact mass: `[m/z] [intensity] "[formula]" [exact_mass]`

The formula matching will utilize the following approach:
1. Apply intensity cutoff to filter peaks
2. For each peak, find the most likely chemical formula that is a subformula of the parent molecule
3. Calculate the exact mass of the assigned formula
4. Account for electron mass and provide correct charges for fragment formulas

This functionality requires both `Molecular_Formula` and `Formula_matching_tolerance` fields to be provided in the input file.

#### Example Output Format with Formula Matching

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
30.033819 2461 "CH4N" 30.034089
55.054611 1497 "C4H7" 55.054770
57.070259 356 "C4H9" 57.070340
68.049652 568 "C5H6" 68.049650
84.080811 2834 "C5H10N" 84.080730
```

Notice how each peak line now includes the assigned formula and its exact mass.

## Current Implementation Limitations

The current version of the tool has the following limitations:

1. **No Formula Matching**: Formula matching for fragments is not yet implemented. This feature is planned for future releases.

2. **Limited Error Recovery**: While the tool logs errors and continues processing other rows, it may not recover gracefully from all potential error conditions.

3. **Basic Metadata Extraction**: Only a subset of available metadata is extracted from mzML files.

4. **Single File Processing**: Each input row creates a separate MSP file; there's no option to combine multiple spectra into a single MSP file.

5. **Limited Mass Accuracy Options**: No support for ppm-based tolerances, only absolute mass tolerances (Da).

## Example

### Sample Input File (TSV):

```
Molecule_name	SMILES	Molecular_formula	Raw_Intensity_Cutoff	Formula_Matching_Tolerance	m/z	RT_seconds	MS_level	mzML_filepath	Spectrum_index
Piperidine	N1CCCCC1	C5H11N	100	0.002	86.09643	305.9	2	/path/to/file.mzML	576
```

### Running the Command:

```bash
pixi run --manifest-path /path/to/pixi.toml python msp_scraper.py input_file.tsv
```

### Output:

The above command will create `Piperidine.msp` in the output directory.

## Future Features

### Formula Matching Implementation

The formula matching functionality will be implemented using optimization techniques to find the best possible formula match for each peak. The key components include:

1. **Formula Decomposition**: Breaking down the parent formula into potential fragment subformulas
2. **Mass Error Minimization**: Finding the formula with the smallest mass error within the specified tolerance
3. **Chemical Plausibility Checks**: Ensuring the assigned formulas make chemical sense
4. **Handling Special Cases**:
   - Accounting for electron mass in charged fragments
   - Supporting isotope patterns
   - Handling adduct ions

Additional dependencies for formula matching will include:
- **molmass**: For handling chemical formulas and isotope calculations
- **cvxpy**: For optimization in formula matching
- **more_itertools**: For efficient generation of subformulas

### Other Planned Enhancements

- **Batch Processing Improvements**: Better handling of large batch processing jobs
- **Extended Metadata Support**: Support for more metadata fields from mzML files
- **Interactive Visualization**: Tools for visualizing extracted spectra
- **Advanced Filtering Options**: More options for filtering and processing spectra
- **Library Building**: Tools for building and searching MSP spectral libraries

## Troubleshooting

### Common Issues

1. **File not found errors**: Ensure all file paths in the input file are correct and accessible.
2. **Index out of range**: Check that the spectrum index exists in the mzML file.
3. **RT not found**: Ensure the retention time is within the range present in the mzML file.
4. **Missing required columns**: Make sure your input file has at least the Molecule_name and mzML_filepath columns.

### Logging

The script generates a log file (default: `msp_scraper.log`) that contains information about the processing steps, warnings, and errors. Check this file for troubleshooting.

## License

MIT

