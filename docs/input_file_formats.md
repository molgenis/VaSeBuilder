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

### Genome reference
As VaSeBuilder can work with both BAM and CRAM alignment files, users need to supply a genome reference in fasta format for each VaSeBuilder run. This genome reference should be the same as the reference used during the mapping process to create the acceptor and donor alignment files.

### Donor FastQ list file
_Note: Using VaSeBuilder-produced FastQ spike-in building blocks is deprecated in favor of using BAMs, as this allows for more flexibility and provides more information._

A donor fastq list file can be provided when running `AssembleValidationSet`. This file requires two columns without a header. The first column should contain the paths to R1 files, the second column the paths to their R2 counter-parts:

```text
/path/to/sample1_R1.fq	/path/to/sample1_R2.fq
/path/to/sample2_R1.fq	/path/to/sample2_R2.fq
/path/to/sample3_R1.fq	/path/to/sample3_R2.fq
```

### Inclusion filter file
A tab-separated inclusion filter file can be provided to indicate which donor variants from which samples to use in the VaSeBuilder run. This is most useful when only certain variants from the donors are required, but when the donor variant files contain many variants.
Variants not specified in this list will therefore be skipped.
The inclusion filter file requires the following 5 columns: sample name/identifier, chromosome name, 1-indexed variant start position, reference allele, and alt allele, with the column headers shown below; additional user-defined custom columns can also be added:

```text
Sample	Chrom	Pos	Ref	Alt	Pathogenicity
sample1	3	12345678	A	T	benign
sample2	4	1357911	GCC	G	likely_pathogenic
sample3	7	24681012	T	TTT	benign
```

### Variant context file
A VaSeBuilder variant context file can serve as the input template file (instead of an acceptor BAM file) when building contexts with `BuildSpikeIns`, and is required for `AssembleValidationSet`. Please see [Variant context file](output_files.md#variant-context-file) in Output Files for its structure.

### Argument file
Argument files can be very useful to run VaSeBuilder. Argument files make it easier to document and redo VaSeBuilder runs. In an argument file, each line should be a command line parameter followed by one or more values (depending on the parameter). Lines starting with a '#' will be ignored. 

```text
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
