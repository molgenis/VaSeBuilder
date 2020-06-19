"""Setuptools configuration file for VaSeBuilder."""
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="VSB_test-tmedina",
    version="0.0.1",
    author="T. Medina",
    author_email="tylerdanmedina@gmail.com",
    description="VSB description.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/molgenis/VaSeBuilder",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: Unix",
        ],
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.18.1",
        "argon2-cffi>=19.2.0",
        "pysam>=0.15.0",
        ]
    )
