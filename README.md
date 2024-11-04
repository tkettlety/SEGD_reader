# SEGD_rev3_reader

A Python package to read seismic data from SEG-D rev 3.0 format.

## Installation:

To install in dev/editable mode (recommended at this stage), run:
```
pip install -e .
```

After, to install optional jupyter lab dependencies (e.g., if using new environment), run:
```
pip install SEGD_rev3_reader[jupyter]
```

## Package contents:

```
SEGD_rev3_reader/
├── SEGD_rev3_reader/
│   ├── __init__.py                    # Initializes package for imports
│   ├── utils.py                       # Utility functions and classes for added functionality
│   ├── segd_reader.py                 # Core reader class and functions
│   └── segd_rev3_csv_headers/                  
│       └── ...                        # CSV files containing header fields from SEG-D rev 3.0 documentation
├── examples/                     
│   ├── sercel_segd3_example.ipynb     # Example notebook for reading Sercel WiNG DFU node data (Costa Rica TAPIR deployment)
│   └── data/                  
│       └── 00102900.segd              # 120 secs Sercel WiNF DFU node data from Costa Rica TAPIR deployment
├── pyproject.toml                     # Build configuration file
├── setup.py                           # Packaging configuration
├── README.md                          # Project description
└── LICENSE                            # License file
```

## License

This package is open-source and licensed under the MIT License.
See the LICENSE file / file headers for more details.

## Version history:

### 0.2 (4th Nov 2024) - Current Version
- Updated reader code to read Stryde node data
- Known issues:
    - Need to identify where GPS location info is stored
    - Does instrument serial number need to be identified
    - Are timings okay? Example file is just over 1 hour.


### 0.1 (1st Nov 2024)
- Initial release of SEGD_rev3_reader
- Only trialled on Sercel WiNG DFU node data
- Limited functionality (reads header/trace data and returns as dictionary or obspy stream)
- Includes additional utility function to rotate 3-component Galperin configuration data to ENZ
