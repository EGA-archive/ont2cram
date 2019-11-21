from setuptools import setup
from ont2cram import META_INFO as _meta

with open("README.md") as f:
    long_description = f.read()

setup(
    name="ont2cram",
    description=_meta['description'],
    long_description=long_description,
    version=_meta['version'],    
    author=_meta['author'],
    install_requires=["h5py", "tqdm", "pysam", "numpy", "parameterized"],
    keywords=_meta['keywords'],
    license=_meta['license'],
    url=_meta['url'],
    classifiers=[        
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3"
    ],    
    py_modules=["ont2cram","cram2ont"],
    scripts=["ont2cram","cram2ont"]
)