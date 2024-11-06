# SEGD_reader

A Python package to read seismic data from SEG-D format.

Currently supports SEG-D revisions 2.1 and 3.0.

## Installation:

To install in dev/editable mode (recommended at this stage), run:
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
│   ├── segd_csv_headers/                  
│       └── ...                          # CSV files containing header fields from SEG-D documentation for different versions
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

## License

This package is open-source and licensed under the MIT License.
See the LICENSE file / file headers for more details.

## Version history:

### 0.3 (6th Nov 2024) - Current Version
- Updated reader code to reader SEG-D revision 2.1 data (tested on SmartSolo node data)
- Renamed package to SEGD_reader to reflect no longer for revision 3.0 only
- Included example notebook for SmartSolo node file

### 0.2 (4th Nov 2024)
- Updated reader code to read Stryde node data
- Possible issues on Stryde data:
    - Need to identify where GPS location info is stored
    - Does instrument serial number need to be identified
    - Are timings okay? Example file is just over 1 hour.


### 0.1 (1st Nov 2024)
- Initial release of SEGD_rev3_reader
- Only trialled on Sercel WiNG DFU node data
- Limited functionality (reads header/trace data and returns as dictionary or obspy stream)
- Includes additional utility function to rotate 3-component Galperin configuration data to ENZ
