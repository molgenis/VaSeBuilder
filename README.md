# VaSeBuilder
Validation Set Builder


## About VaSeBuilder
### Short introduction
Often when a computational analysis consisting of multiple steps is run routinely, a pipeline is constructed. The pipeline connects the individual analysis steps and automates the execution of each and transition to the next. Pipelines can be evaluated to determine whether each step and transition to the next works correctly. Simultaneously, validation of the pipeline to determine whether it can indeed provide relevant results is equally important. Pipelines can be validated in separate and isolated runs but this requires time. We aim to combine the running of analyses and validation of the pipeline.\
VaSeBuilder (Validation Set Builder) can be used to construct two FastQ files with forward and reverse reads for the validation of the Molgenis NGS_DNA pipeline. To construct these two FastQ files data from samples already processed with the NGS_DNA pipeline are used to modify a template.\
The sample data should consist of a BAM file (containing aligned reads) and a VCF file containing identified variants.\
The template can for example be the NA12878 sample and should first be processed with the NGS_DNA pipeline. VaSeBuilder only requires the two FastQ and the produced BAM file.\
\


### What does VaSeBuilder do?
For each provided sample, VaSeBuilder extracts BAM reads overlapping with a variant noted in the VCF file. The mate of the overlapping read is also included. From these reads leftmost and rightmost positions are determined. These two positions constitute the context start and stop. Variants located within a previously established context are skipped. From the template BAM, reads (including their mates) overlapping with the variant are also extracted.\
Once all samples have been processed, the two template FastQ files are processed to produce the two new validation FastQ files. Template reads not overlapping with any variant are written to the validation FastQ files. Template reads overlapping with a variant are replaced with sample reads overlapping the same variant.\
This produces two FastQ files for which is know which variants they contain and therefore which variants the pipeline should be able to identify.\
Currently (feb. 2019) VaSeBuilder only works with 'simple' genomic variants such as SNPs and small indels, but this may very well be expanded in the future :)\
\



## Basic usage
To run VaSeBuilder run vase.py and set all required parameters via the command line. The list of parameters is listed below and can also be viewed via _python vase.py -h_. VaSeBuilder requires the location of one or more valid VCF and BAM directories. If one of the directories does not exist of does not contain any VCF files the folder will be skipped but the program can still continue if data from other directories can still be used.\
When running, the program will output information about the actions it is performing and potential problems it may have encountered.


### Program parameters
#### Required parameters
* __-v__/__--donorvcf__: Provide the path to one or more valid directory containing VCF files. **For example:** *--vcfin /vcfData/vcfDirectory1 /vcfData/vcfDirectory2*
* __-b__/__--donorbam__: Provide the path to one or more valid directory containing BAM files. **For example:** *--bam /bamData/bamDirectory1 /bamData/bamDirectory2*
* __-a__/__--acceptorbam__: Provide the location of the BAM file of the sample that will be used as the template to create the validation FastQ files from."**For example:** *--templatebam /templateData/template.bam*
* __-1__/__--templatefq1__: Provide the location of the first FastQ file that will be used as the template to produce the first validation FastQ file. **For example:** *--templatefq1 /fqData/template_reads_1.fastq.gz*
* __-2__/__--templatefq2__: Provide the location of the first FastQ file that will be used as the template to produce the second validation FastQ file. **For example:** *--templatefq2 /fqData/template_reads2.fastq.gz*
* __-o__/__--out__: Provide the location where to write the output files to.

#### Optional parameters
* __-of__/__--fastqout__: Provide the name VaSeBuilder should use a name prefix for the FastQ files. **For example:** *--fastqout /outData/VaSeFq*
* __-ov__/__--varcon__: Provide the file name VaSeBuilder should write the used variant contexts to. **For example:** *--varcon /outData/variant_contexts.txt*
* __-od__/__--donorbread__: Provide the name and location where VaSeBuilder should write BAM reads associated to VCF variants **For example:** *--varbead /outData/variants_bamreads.txt*
* __-oa__/__--acceptorbread__: Provide the name and location where VaSeBuilder should write template BAM reads associated with VCF variats to **For example:** *--templatebread /outData/template_variant_bamreads.txt*
* __-l__/__--log__: You can provide the name and location where to write the log file to. This parameter will write the log file to the current working directory if not set. **For example:** *--log /outData/vaselog.log*


### Program output
Aside from the gzipped forward and reverse read FastQ files, the program also outputs a set of text files. These files are all related to the variant contexts. Although the names of most files can be set via the optional parameters, below we used the default names to describe each
* __acceptorbread.txt__: Tab separated file containing variant 
* __donorbams.txt__: Contains the list of BAM files as a read donor per sample separated by tabs. The completed paths to each BAM file are reported and separated from each other with ' ; ' (including the spaces).
* __donorbread.txt__: Contains the list of donor BAM reads per variant context. Donor reads are separated from each other with ' ; ' (including the spaces)
* __donorvcfs.txt__: Contains the list of VCF files as a read donor per sample separated by tabs. The completed paths to each VCF file are reported and separated from each other with ' ; ' (including the spaces).
* __unmappedmatereads.txt__: Contains the read identifiers of reads with an unmapped mate for each sample. Acceptor reads with unmapped mates are also 
* __varcon.txt__: Contains the essential data, consisting of VCF variant Id, Sample, context chromosome, context start position and context stop position, for each identified and processed variant context separated by tabs.



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



## Still in progress
(As of feb. 14th, 2019)
* __Making proper unittests:__ _Although there are unittests, these should still be improved and a few still need to be added._
* __Genuine first functional test:__ _Does the two newly created fastq files make sense (are the reads replaced correctly, what about coverage, what about QC differences after replacing reads, etc)._



## Future things
VaSeBuilder might change into a python project installable via pip.
VaSeBuilder might very well be updated and changed to work with more complex variants
Make VaSeBuilder into a multiprocessor program to speed up te process.
Perhaps add the number of donor and acceptorr BAM reads as columns of the varcon.txt file.

Please let me know if this documentation is missing things or which things are unclear.