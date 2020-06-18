"""Setuptools configuration file for VaSeBuilder."""
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="example_package-tmedina",
    version="0.0.1",
    author="T. Medina",
    author_email="tylerdanmedina@gmail.com",
    description="Test description.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/molgenis/VaSeBuilder",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3",
        "Operating System :: Unix",
        ],
    python_requires=">=3.7",
    )
