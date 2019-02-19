# VaSeBuilder
Validation Set Builder


## About VaSeBuilder
### Short introduction
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
To run VaSeBuilder run vase.py and set all required parameters via the command line. The list of parameters is listed below and can also be viewed via _python vase.py_. VaSeBuilder requires the location of one or more valid VCF and BAM directories. If one of the directories does not exist of does not contain any VCF files the folder will be skipped but the program can still continue if data from other directories can still be used.\
When running, the program will output information about the actions that it is performing and potential problems it encountered. This program also 



### Program parameters
* __--vcfin__: Provide the path to one or more valid directory containing VCF files. **For example:** *--vcfin /vcfData/vcfDirectory1 /vcfData/vcfDirectory2*
* __--bamin__: Provide the path to one or more valid directory containing BAM files. **For example:** *--bam /bamData/bamDirectory1 /bamData/bamDirectory2*
* __--templatebam__: Provide the location of the BAM file of the sample that will be used as the template to create the validation FastQ files from."**For example:** *--templatebam /templateData/template.bam*
* __--templatefq1__: Provide the location of the first FastQ file that will be used as the template to produce the first validation FastQ file. **For example:** *--templatefq1 /fqData/template_reads_1.fastq.gz*
* __--templatefq2__: Provide the location of the first FastQ file that will be used as the template to produce the second validation FastQ file. **For example:** *--templatefq2 /fqData/template_reads2.fastq.gz*
* __--fastqout__: Provide the directory where VaSeBuilder should write the two new FastQ files to. **For example:** *--fastqout /outData*
* __--varcon__: Provide the name and location where VaSeBuilder should write the used variants and their contexts to. **For example:** *--varcon /outData/variants_contexts.txt*
* __--varbread__: Provide the name and location where VaSeBuilder should write BAM reads associated to VCF variants **For example:** *--varbead /outData/variants_bamreads.txt*
* __--templatebread__: Provide the name and location where VaSeBuilder should write template BAM reads associated with VCF variats to **For example:** *--templatebread /outData/template_variant_bamreads.txt*
* __--log__: You can provide the name and location where to write the log file to. This parameter is not required and will write the log file to the current working directory if not set. **For example:** *--log /outData/vaselog.log*



## Questions & Answers
**Q: Which python version is required?**\
**A:** VaSeBuilder has been written, run and tested in Python 3. To run the program on Python 2 the code might need to be changed first.

**Q: Which python libraries are required?**\
**A:** VaSeBuilder uses _pysam_, _gzip_, _Biopython_ and _numpy_. The easiest way to obtain these libraries is to install them via pip (pip install _library_). Note that the pysam library itself requires htslib which doesn't run on the Windows platform.

**Q: Can VaSeBuilder be run on Windows?**\
**A:** VaSeBuilder might run on Windows by replacing the pysam library with the bamnostic library that can also be used to work with BAM files but does not require htslib.

**Q: Are sample identifiers required?**\
**A:** Yes, VaSeBuilder uses the sample identifiers in BAM and VCF files to link the two files. You can check if a BAM has sample information by examing the header for a line similar to:
_@RG	ID:Sample	SM:Sample_
You can extract the header (samtools view -H bamFile.bam > bamHeader.txt), add a sample by adding the above line with the correct sample name and then reheader the BAM file (samtools reheader bamHeader.txt bamFile.bam > BamFileWithSample.bam).


## Still in progress
(As of feb. 14th, 2019)
* __Making proper unittests:__ _Although there are unittests, these should still be improved and a few still need to be added._
* __Genuine first functional test:__ _Does the two newly created fastq files make sense (are the reads replaced correctly, what about coverage, what about QC differences after replacing reads, etc)._


## Future things
VaSeBuilder might change into a python project installable via pip.
VaSeBuilder might very well be updated and changed to work with more complex variants
