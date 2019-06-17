from setuptools import setup

with open("README.md") as f:
    long_description = f.read()

setup(
    name="ont2cram",
    description="Oxford Nanopore HDF/Fast5 to CRAM conversion tool",
    long_description=long_description,
    version = "0.0.1",
    author="EMBL-EBI",
    install_requires=["h5py", "tqdm", "pysam", "numpy", "parameterized"],
    keywords=["ONT", "FAST5", "HDF", "CRAM"],
    license="Apache License, Version 2.0",
    url="https://github.com/EGA-archive/ont2cram",
    classifiers=[        
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3"
    ],    
    py_modules=["ont2cram","cram2ont"],
    scripts=["ont2cram","cram2ont"]
)