# SEGD_reader

A Python package to read seismic data from SEG-D format.

Currently supports SEG-D revisions 2.1 and 3.0.

## Installation:

To install in dev/editable mode (recommended at this stage), `cd` to the repository on your local storage and run:
```
pip install -e .
```

After, to install optional jupyter lab dependencies (e.g., if using new environment), run:
```
pip install SEGD_reader[jupyter]
```

## Package contents:

```
SEGD_reader/
├── SEGD_reader/
│   ├── __init__.py                      # Initializes package for imports
│   ├── utils.py                         # Utility functions and classes for added functionality
│   ├── segd_reader.py                   # Core reader class and functions
│   ├── json/                  
│       └── ...                          # JSON files containing header fields from SEG-D documentation for different versions
├── examples/                     
│   ├── sercel_rev3_0_example.ipynb      # Example notebook for reading Sercel WiNG DFU node data (Costa Rica TAPIR deployment)
│   ├── stryde_rev3_0_example.ipynb      # Example notebook for reading Stryde node data (which deployment?)
│   ├── smartsolo_rev2_1_example.ipynb   # Example notebook for reading SmartSolo node data (which deployment?)
│   └── data/                  
│       └── ...                          # Example SEG-D files
├── pyproject.toml                       # Build configuration file
├── setup.py                             # Packaging configuration
├── README.md                            # Project description
└── LICENSE                              # License file
```

## Known issues:
- Sensitivity, lat, lon, etc stored in manufacturer-defined trace extension headers in SEG-D 2.1 files. Only know this structure for Sercel SEG-D files.
- Need to find polarity correction too (sometimes positive sometimes negative)
- License needs updating to reflect non-commerical use and to credit initial source of revision 2.1 csv files

## License

This package is open-source and licensed under the MIT License.
See the LICENSE file / file headers for more details.

## Version history:

### 0.3.4 (19th Jun 2025) - Current Version
- Extracts more header info from Sercel SEG-D 2.1 files (sensitivity, lat, lon, serial number, ...)
- Trace stats in obspy more explicitly say whether descale multiplier and sensitivity corrections have been applied

### 0.3.3 (19th Jun 2025)
- Replaced header spec csv files with a single JSON file for each SEG-D version and updated reader class to use these
- Made SEG_D_to_stream a little more stable with respect to handling empty traces and converting traces to int

### 0.3.2 (24th Jan 2025)
- Updated SEG_D_to_stream function to apply descale multiplier from trace header
- Updated checks when trying to rotate from Galperin configuration to ENZ (utils)

### 0.3.1 (8th Nov 2024)
- Fixed reader code to read all traces in SEG-D revision 2.1 data
- Updated SmartSolo example notebook to reflect above changes
- Tested on SmartSolo and Sercel SEG-D revision 2.1 data

### 0.3.0 (6th Nov 2024)
- Updated reader code to read SEG-D revision 2.1 data (tested on SmartSolo node data)
- Renamed package to SEGD_reader to reflect no longer for revision 3.0 only
- Included example notebook for SmartSolo node file

### 0.2.0 (4th Nov 2024)
- Updated reader code to read Stryde node data
- Possible issues on Stryde data:
    - Need to identify where GPS location info is stored
    - Does instrument serial number need to be identified
    - Are timings okay? Example file is just over 1 hour.


### 0.1.0 (1st Nov 2024)
- Initial release of SEGD_rev3_reader
- Only trialled on Sercel WiNG DFU node data
- Limited functionality (reads header/trace data and returns as dictionary or obspy stream)
- Includes additional utility function to rotate 3-component Galperin configuration data to ENZ
