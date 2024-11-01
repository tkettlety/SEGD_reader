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

SEGD_rev3_reader/
├── SEGD_rev3_reader/
│   ├── __init__.py                    # Initializes package for imports
│   ├── utils.py                       # Utility functions and classes for added functionality
│   ├── segd_reader.py                 # Core reader class and functions
│   └── segd_rev3_csv_headers/                  
│       └── ...                        # CSV files containing header fields from SEG-D rev 3.0 documentation
├── examples/                     
│   └── sercel_segd3_example.ipynb     # Example notebook for reading Sercel WiNG DFU node data (Costa Rica TAPIR deployment)
├── pyproject.toml                     # Build configuration file
├── setup.py                           # Packaging configuration
├── README.md                          # Project description
└── LICENSE                            # License file