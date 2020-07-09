"""Setuptools configuration file for VaSeBuilder."""
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="VaSeBuilder",
    version="20.7.1",
    author="T. Medina",
    author_email="tylerdanmedina@gmail.com",
    description="VaSeBuilder, a bioinformatic tool for artificially combining multiple NGS samples into a single, hybrid sample.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/molgenis/VaSeBuilder",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: Unix",
        ],
    python_requires=">=3.6",
    install_requires=[
        "numpy>=1.18.1",
        "argon2-cffi>=19.2.0",
        "pysam>=0.15.0",
        ],
    # scripts=["bin/VaSeBuilder"],
    entry_points={
        "console_scripts": ["VaSeBuilder = VaSeBuilder.__main__:main"]
        }
    )
