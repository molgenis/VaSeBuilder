# VaSeBuilder input file formats

## About input files
To run VaSeBuilder specific VaSeBuilder input files might need to be provided depending on the selected run mode.


### Donor variant and donor alignment list files


### Donor fastq list file
A donor fastq list file needs to be provided when running VaSeBuilder AC-mode. Currently (September 23rd, 2019) this 
file requires two columns without a header for. The first column should contain the paths to R1 files, the second column
 the paths to R2 files.


sample1_R1.fq &emsp;sample1_R2.fq<br />
sample2_R1.fq &emsp;sample2_R2.fq


### Variant list file


### Variant context file
A VaSeBuilder variant context file can serve as an input file when running AC, DC, FC or PC mode. Please see the output 
files description for the variant context file.


### Configuration file
Configuration files allow users to save the VaSeBuilder run parameters and values in a text file. This makes rerunning 
analyses on the same data, as well as sharing the run parameters much easier. Parameters and values can be specified as 
PARAMETERNAME=value. Parameter names do not need to be capitalized, entries such as parametername = value are also 
valid. Parameters templatefq1 and templatefq2 can have multiple values separated by a comma: 
templatefq1 = testdata/fastqs/acceptor1_1.fq , testdata/fastqs. (The space around the comma is optional)
Comments explaining or describing the use of the configuration file can be added by starting a line with #.

Valid configuration file parameter names:
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


To run VaSeBuilder with a configuration file use: python vase.py –c path/to/con fig/file.cfg
Please see the example config files in VaSeBuilder/config/examples for examples per run mode
