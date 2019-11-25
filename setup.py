from setuptools import setup
import ont2cram as package

with open("README.md") as f:
    long_description = f.read()

setup(
    name=package.__name__,
    description=package.__description__,
    long_description = long_description,
    long_description_content_type="text/markdown",
    version=package.__version__,
    author=package.__author__,
    install_requires=[
        'h5py>=2.10',
        "ont_fast5_api>=2.0.0",
        "tqdm>=4.39.0",
        "pysam>=0.15.3",
        "numpy>=1.17.4 ",
        "parameterized>=0.7.1"],
    keywords=package.__keywords__,
    license=package.__license__,
    url=package.__url__,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3"
    ],
    package_dir = {package.__name__: package.__name__},
    package_data = {package.__name__: ['test/test_data/*', 'test/test_data/*/*']},
    entry_points = {'console_scripts': ['ont2cram=ont2cram.__main__:main']}
)
