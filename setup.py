from setuptools import setup, find_packages

setup(
    name="SEGD_reader",
    version="0.3.2",
    description="A Python package to read seismic data from SEG-D format",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/SEGDreader/segd_reader",
    license="MIT",
    packages=find_packages(exclude=["examples", "extra", "build", "dist", "segd_docs"]),
    package_data={
        "SEGD_reader": ["segd_csv_headers/**/*"],
    },
    include_package_data=True,
    install_requires=[
        "numpy",
        "obspy"
    ],
    extras_require={
        "jupyter": ["jupyterlab", "ipython"]
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        'Topic :: Scientific/Engineering :: Physics'
    ],
    python_requires='>=3.7',
)