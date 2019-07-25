# VaSeBuilder
Validation Set Builder
&nbsp;

&nbsp;

## About VaSeBuilder
### Short introduction
Often when a computational analysis consisting of multiple steps is run routinely, a pipeline is constructed. The pipeline connects the individual analysis steps and automates the execution of each and transition to the next. Pipelines can be evaluated to determine whether each step and transition to the next works correctly. Simultaneously, validation of the pipeline to determine whether it can indeed provide relevant results is equally important. Pipelines can be validated in separate and isolated runs but this requires time. We aim to combine the running of analyses and validation of the pipeline.\
VaSeBuilder (Validation Set Builder) can be used to construct two FastQ files with forward and reverse reads for the validation of the Molgenis NGS_DNA pipeline. To construct these two FastQ files data from samples already processed with the NGS_DNA pipeline are used to modify a template.\
The sample data should consist of a BAM file (containing aligned reads) and a VCF file containing identified variants.\
The template can for example be the NA12878 sample and should first be processed with the NGS_DNA pipeline. VaSeBuilder only requires the two FastQ and the produced BAM file.

### What does VaSeBuilder do?
For each sample, VaSeBuilder iterates over the variants within the VCF/BCf file (if a variant list is provided only 
variants satisfying that list will be used). First reads overlapping directly with the variant are identified in both 
the acceptor/template and the donor BAM/CRAM file of the same sample. From these reads and their mates the left and 
right most position is determined which establishes the acceptor context and donor context respectively. From these two 
contexts the variant context is determined. This context spans the absolute left and right most genomic positions of 
both contexts and can (quite often) be larger than either. Reads overlapping with this variant context and teir mates 
are then identified and saved.\
After processing all donor samples, the acceptor/template fastq files are used to produce the validation set fastq 
files. Acceptor reads within a variant contexts are excluded from these new fastq files and are replaced with donor 
reads that are located in the same variant context. This produces a set of fastq files for with known variants that 
should be able to be identified and which reads carry those variants.\
Currently VaSeBuilder only works with 'simple' genomic variants such as SNPs and small indels, but this may very well 
be expanded in the future :)


### What are acceptor, donor and variant contexts?
VaSeBuilder creates the validation set by identifying for each variant which acceptor reads to be exchanged with which 
donor reads. The acceptor context is the window established by the leftmost and rightmost genomic positions of reads 
directly overlapping with the variants and their read mates (which likely do not overlap). The donor context works the 
same as the acceptor context but instead uses a donor file.\
The variant context is established by combining the minimum and maximum border of acceptor and donor context and thereby
 spans both. Since the variant context spans a larger area than both acceptor and donor context individually, acceptor 
and donor reads and their mates are identified again. Acceptor reads overlapping with the variant context and their 
mates are excluded and replaced by donor reads overlapping with the variant context.


### Required software
* Python 3.6 or higher
* Pysam 0.14 or higher
* Linux file command v5.37 or higher
* HTSlib 1.7 or higher
* SAMtools 1.7 or higher


### Important to know
VaSeBuilder is intended to build a validation set from acceptor and donor data that was sequenced with the same
sequencer/sequencing platform and treated with the same preparation and capturing kit. (Please see the documentation
 later as to why...)
&nbsp;

&nbsp;

## Basic usage
To run VaSeBuilder run vase.py and set all required parameters via the command line. The program parameters are listed 
below and can also be viewed via _python vase.py -h_. VaSeBuilder requires the location of one or more valid VCF and BAM
 directories. If one of the directories does not exist of does not contain any VCF files the folder will be skipped but 
the program can still continue if data from other directories can still be used.\
When running, the program will output information about the actions it is performing and potential problems it may have encountered.\
__To run the program for example:__
python vase.py -v /data/vcf_file_list.txt -b /data/bamcram_file_list.txt -a /data/acceptor_sample.bam 
-1 /data/acceptor_R1.fq.gz -2 /data/acceptor_R2.fq.gz -o /data/output -r /data/human_reference.fa


## Run modes
VaSeBuilder offers several run mode. By default, VaSeBuilder will run the 'normal' way as described above. In total 
VaSeBuilder offers five different running modes that can be specified using the -m/--runmode parameter.
* __A__: Existing donor fastq files specified in a list file are added to acceptor/template fastq files.
* __D__: Variant contexts are identified and written to a variant context file. Donor read data for all variant contexts
 is written to fastq files.
* __F__: Variant contexts are identified and written to a variant context file. Acceptor reads are filtered from the 
acceptor/template fastq files and donor reads are added.
* __P__: Variant contexts are identified and written to a variant context file. Donor read data for each variant context
 is written to separate fastq files.
* __X__: Variant contexts with acceptor and donor reads are identified and written to a variant context file but no 
fastq files are outputted.
Additionally one of the five runmodes, a 'C' can be added to continue from an existing variant context file. If the 
additional 'C' is used, the -vc/--varconin parameter is required.

### Required parameters per run mode
* __A__: -1/--templatefq1 , -2/--templatefq2 , -o/--out , -dq/--donorfastqs , -vc/--varconin
* __D__: -v/--donorvcf , -b/--donorbam , -a/--acceptorbam , -r/--reference , -o/--out
* __F__: -v/--donorvcf , -b/--donorbam , -a/--acceptorbam , -r/--reference , -o/--out , -1/--templatefq1 , 
-2/--templatefq2 , 
* __P__: -v/--donorvcf , -b/--donorbam , -a/--acceptorbam , -r/--reference , -o/--out
* __X__: -v/--donorvcf , -b/--donorbam , -a/--acceptorbam , -r/--reference , -o/--out


### Program parameters
#### Required run mode parameters
* __-v__/__--donorvcf__: File containing the locations of donor VCF/BCF files to use (one VCF/BCF file per line). For example: *-v /data/vcf_list_file.txt*
* __-b__/__--donorbam__: File containing the locations of donor BAM/CRAM files to use (one BAM/CRAM file per line). **For example:** *--donorbam /bamData/bamDirectory1 /bamData/bamDirectory2*
* __-a__/__--acceptorbam__: BAM file of the sample that will be used as the template to create the validation FastQ files from."**For example:** *--acceptorbam /templateData/template.bam*
* __-1__/__--templatefq1__: Provide the location of the first FastQ file that will be used as the template to produce the first validation FastQ file. For example: *-1 /fqData/template_reads_R1.fastq.gz*
* __-2__/__--templatefq2__: Provide the location of the first FastQ file that will be used as the template to produce the second validation FastQ file. For example: *-2 /fqData/template_reads_R2.fastq.gz*
* __-o__/__--out__: Provide the location where to write the output files to.
* __-r__/__--reference__: Provide the reference genome in fasta format used for mapping (This reference willk be used for all BAM/CRAM files). For example: *-r /data/human_reference.fa*
* __-dq__/__--donorfastqs__: Provide a file containing locations to donor fastq files to use for constructing the validation set. This list file should contain two columns; the first containing the R1 fastq files, the second the R2 fastq files

#### Optional parameters
* __-of__/__--fastqout__: Provide the name VaSeBuilder should use a name prefix for the FastQ files. **For example:** *--fastqout /outData/VaSeFq*
* __-ov__/__--varcon__: Provide the file name VaSeBuilder should write the used variant contexts to. **For example:** *--varcon /outData/variant_contexts.txt*
* __-l__/__--log__: You can provide the name and location where to write the log file to. This parameter will write the log file to the current working directory if not set. **For example:** *--log /outData/vaselog.log*
* __-!__/__--debug__: Run the program in debug mode. (This will create a more detailed log and additional output files)
* __-vl__/__--variantlist__: List of specific variants to use to build the new Validation/Variant set

### Program output
Aside from the forward and reverse read FastQ files, the program also outputs a set of text files. These files are all 
related to the variant contexts. Although the names of some files can be set via the optional parameters, below we used 
the default names to describe each.

#### Default output files
* __donorbams.txt__: List of donor BAM/CRAM files used to create the validation set
* __donorvcfs.txt__: List of donor VCF/BCF files used to create the validation set
* __varcon.txt__: List of variant contexts that build the validation set
* __varconstats.txt__: Basic variant context statistics (average and median)
* __variantcontexts.bed__: The variants contexts in bed file format.
* __VaSeBuilder.log__: The log file containing data about the run of the program

#### Debug output files
* __acceptorcontexts.txt__: List of established acceptor contexts.
* __acceptorcontextstats.txt__: Basic acceptor context statistics
* __acceptor_unmapped.txt__: List of acceptor read identifier per acceptor context that have unmapped read mates 
* __acceptor_positions.txt__: List of all acceptor R1 read left starting positions and R2 read right ending positions per acceptor context
* __donorcontexts.txt__: List of established donor contexts.
* __donorcontextstats.txt__: Basic donor context statistics 
* __donor_unmapped.txt__: List of donor read identifiers per donor context that have unmapped read mates.
* __donor_positions.txt__: List of all donor R1 read left starting positions and R2 read right ending positions per donor context
* __varcon_unmapped_acceptor.txt__: List of variant context acceptor read identifiers that have unmapped mates per context
* __varcon_unmapped_donor.txt__: List of variant context donor read identifiers that have unmapped mates per context
* __varcon_positions_acceptor.txt__: List of all variant context acceptor R1 read left starting positions and R2 read right ending positions per variant context
* __varcon_positions_donor.txt__: List of all variant context donor R1 read left starting positions and R2 read right ending positions per variant context
&nbsp;

&nbsp;


## Questions & Answers
### General questions
**Q: Which python version is required?**\
**A:** VaSeBuilder has been created to run in Python 3. To run the program on Python 2 the code might need to be changed first.

**Q: Which python libraries are required?**\
**A:** VaSeBuilder uses _pysam_ and _gzip_. The easiest way to obtain these libraries is to install them via pip (pip install _library_). Note that the pysam library itself requires htslib which doesn't run on the Windows platform.

**Q: Which pysam version is required?**\
**A:** Pysam 0.14.0 or higher is required as VaSeBuilder uses the get_forward_sequence() and get_forward_qualities() methods that were introduced in pysam 0.14.0.

**Q: Can VaSeBuilder be run on Windows?**\
**A:** Some adjustments would need to be made as pysam can't be used on Windows. The python library bamnostic can replace pysam for working with BAM files. Reading and processing VCF would need to be done with the default python file reading/writing functionality or other libraries (if available).

**Q: Are sample identifiers required?**\
**A:** Yes, VaSeBuilder uses the sample identifiers in BAM and VCF files to link the two files. You can check if a BAM has sample information by examing the header for a line similar to:
_@RG	ID:Sample	SM:Sample_
You can extract the header (samtools view -H bamFile.bam > bamHeader.txt), add a sample by adding the above line with the correct sample name and then reheader the BAM file (samtools reheader bamHeader.txt bamFile.bam > BamFileWithSample.bam).


### Specific questions
**Q: Is there a simple way to check whether the donor and acceptor reads are indeed added and removed respectively?**\
**A:** Currently I'm working on a set of scripts, as well as others, that I have called VaSeUtils which can be used to check this.

**Q: What is the effect of using data from different sequencers, prep-kits and/or capturing kits?**\
**A:** To still be investigated...
&nbsp;

&nbsp;


## Future things
VaSeBuilder might change into a python project installable via pip.\
VaSeBuilder might very well be updated and changed to work with more complex variants.\
Make VaSeBuilder into a multiprocessor program to speed up te process.\
Let VaSeBuilder work with more complex variants.\
Let VaSeBuilder work with trio data.

Please let me know if this documentation is missing things or which things are unclear.