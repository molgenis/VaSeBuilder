# VaSeBuilder input file formats

## About input files
To run VaSeBuilder, specific VaSeBuilder input files might need to be provided depending on the selected run mode.<br />
<br />


### Donor variant and donor alignment list files
The donor variant list file should contain a list of paths to a variant, VCF or BCF, file with one path per file line.
The donor alignment list file should be similar and contain one path to a donor alignment file per file line.
The files listed will be used as donor variant and alignments files respectively when constructing variant contexts.
<br /><br />


### Donor fastq list file
A donor fastq list file needs to be provided when running VaSeBuilder AC-mode. Currently (September 23rd, 2019) this 
file requires two columns without a header. The first column should contain the paths to R1 files, the second column
 the paths to R2 files:

_/path/to/sample1_R1.fq &emsp;/path/to/sample1_R2.fq<br />
/path/to/sample2_R1.fq &emsp;/path/to/sample2_R2.fq<br />
/path/to/sample3_R1.fq &emsp;/path/to/sample3_R2.fq_

Each row in this file therefore represents a single sample.<br /><br />


### Variant list file
A variant list file can be provided as a filter to indicate which donor variants to use in the VaSeBuilder run. This is 
most useful when only certain variants from the donors are required but the donor variant files contain many variants.
Variants not specified in the variant list will therefore be skipped.
The variant list file needs at least three columns: sample name/identifier, chromosome name, variant position:

_Sample1 &emsp;1 &emsp;101<br />
Sample2 &emsp;2 &emsp;202<br />_

(In the future two more columns might be required: variant reference allele, variant alternative allele)<br /><br />


### Variant context file
A VaSeBuilder variant context file can serve as an input file when running AC, DC, FC or PC mode. Please see the output 
files description for the variant context file for it's structure.<br /><br />


### Configuration file
Configuration files allow users to save the VaSeBuilder run parameters and values in a text file. This makes rerunning 
analyses on the same data, as well as sharing the run parameters much easier. Parameters and values can be specified as 
PARAMETERNAME=value. Parameter names do not need to be capitalized, entries such as parametername = value are also 
valid. Parameters templatefq1 and templatefq2 can have multiple values separated by a comma: 
templatefq1 = testdata/fastqs/acceptor1_1.fq , testdata/fastqs/acceptor1_2.fq. (The space around the comma is optional)
Comments explaining or describing the use of the configuration file can be added by starting a line with #.

<br /><u>Valid configuration file parameter names:</u>
* __RUNMODE:__ The run mode of VaSeBuilder to use 
* __DONORVCF:__ Path to file with list of donor variant files to use
* __DONORBAM:__ Path to file with list of donor alignment files to use
* __ACCEPTORBAM:__ Path the alignment file to use as acceptor
* __TEMPLATEFQ1:__ Path(s) to R1 fastq file(s) to use as template for building the validation set
* __TEMPLATEFQ2:__ Path(s) to R2 fastq file(s) to use as template for building the validation set
* __REFERENCE:__ Path to the fasta genome reference
* __OUT:__ Path to directory to write output files to
* __LOG:__ Path to write log file to
* __DEBUG:__ Run debug mode (valid values: T, True, true, TRUE, 1).
* __VARIANTLIST:__ Path to file containing  variants to use
* __VARCONIN:__ Path to a variant context file
* __FASTQOUT:__ Name prefix for validation fastq files
* __VARCON:__ Name prefix for variant context output file

To run VaSeBuilder with a configuration file use: python vase.py –c path/to/config/file.cfg
Please see the example config files in VaSeBuilder/config/examples for examples per run mode
