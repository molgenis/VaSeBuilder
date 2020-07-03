# VaSeBuilder input file formats

## About input files
To run VaSeBuilder, specific VaSeBuilder input files might need to be provided depending on the selected run mode.<br />
<br />


### Donor variant and donor alignment list files
The donor variant list file should contain a list of paths to a variant, VCF or BCF, file with one path per file line.
The donor alignment list file should be similar and contain one path to a donor alignment file per file line.
The files listed will be used as donor variant and alignments files respectively when constructing variant contexts.
Only paths to the variant/alignment files should be provided, not the index files. These are assumed to be in the same 
directory with the same name, i.e.: sample.bam & sample.bam.bai
<br /><br />

### Genome reference
As VaSeBuilder can work with both BAM and CRAM alignment files, users need to supply a genome reference in fasta format for each VaSeBuilder run. This genome reference should be the same as the reference used during the mapping process to create the acceptor and donor alignment files.


### Donor FastQ list file
A donor fastq list file needs to be provided when running VaSeBuilder AC-mode. Currently (September 23rd, 2019) this 
file requires two columns without a header. The first column should contain the paths to R1 files, the second column
 the paths to R2 files:

```
/path/to/sample1_R1.fq &emsp;/path/to/sample1_R2.fq
/path/to/sample2_R1.fq &emsp;/path/to/sample2_R2.fq
/path/to/sample3_R1.fq &emsp;/path/to/sample3_R2.fq
```

Each row in this file therefore represents a single sample.
<br /><br />


### Variant list file
A variant list file can be provided as a filter to indicate which donor variants to use in the VaSeBuilder run. This is 
most useful when only certain variants from the donors are required but the donor variant files contain many variants.
Variants not specified in the variant list will therefore be skipped.
The variant list file needs at least three columns: sample name/identifier, chromosome name, variant position:

```
Sample1 &emsp;1 &emsp;101  
Sample2 &emsp;2 &emsp;202  
```
<br />

(In the future two more columns might be required: variant reference allele, variant alternative allele)
<br /><br />


### Variant context file
A VaSeBuilder variant context file can serve as an input file when running AC, DC, FC or PC mode. Please see the output 
files description for the variant context file for it's structure.
<br /><br />


### Argument file
Argument files can be very useful to run VaSeBuilder. Argument files make it easier to document and redo VaSeBuilder runs. In an argument file, each line should be a command line parameter followed by one or more values (depending on the parameter). Lines starting with a '#' will be ignored. 
<br />
```
BuildValidationSet  
--acceptor-bam GIAB_NA12878.bam  
--donor-vcf-list MyVCFList.txt  
--donor-bam-list MyBAMList.txt  
--reference hg37.fa  
--inclusion-filter MyFilterList.txt  
--acceptor-fq-r1 GIAB_L1_1.fq.gz GIAB_L2_1.fq.gz  
--acceptor-fq-r2 GIAB_L1_2.fq.gz GIAB_L2_2.fq.gz  
--fastq-out MyVaseSet  
--log MyLog.log  
#--no-hash
```
