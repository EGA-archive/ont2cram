from setuptools import setup
import ont2cram as pkg

with open("README.md") as f:
    long_description = f.read()

setup(
    name=pkg.__name__,
    description=pkg.__description__,
    long_description = long_description,
    long_description_content_type="text/markdown",
    version=pkg.__version__,
    author=pkg.__author__,
    install_requires=[
        'h5py>=2.10',
        "ont_fast5_api>=2.0.0",
        "tqdm>=4.39.0",
        "pysam>=0.15.3",
        "numpy>=1.17.4 ",
        "parameterized>=0.7.1"],
    keywords=pkg.__keywords__,
    license=pkg.__license__,
    url=pkg.__url__,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3"
    ],
    packages = [pkg.__name__],
    package_dir = {pkg.__name__: pkg.__name__},
    package_data = {pkg.__name__: ['test_data/*', 'test_data/*/*']},
    entry_points = {'console_scripts': ['ont2cram=ont2cram.__main__:main']}
)
