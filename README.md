# VaSeBuilder
Validation Set Builder


## About VaSeBuilder
VaSeBuilder can be used to construct two FastQ files with forward and reverse reads that can be used for validation of the Molgenis NGS_DNA pipeline. 

## Basic usage
To run VaSeBuilder run vase.py and set all required parameters. VaSeBuilder requires the location of one or more valid VCF

### Program parameters
--vcfin: Provide at one or more valid folder(s) containing VCF files. For example: --vcfin /vcfData/vcfFolder1 /vcfData/vcfFolder2
--bamin: 
--valbam: 

##Future things
VaSeBuilder might change into a python project installable with pip.
