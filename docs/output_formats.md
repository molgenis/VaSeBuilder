# VaSeBuilder output file formats

## Standard output files
VaSeBuilder runs output several output files containing the essential data about the run. Some files are only outputted 
when running specific modes.

### Donor alignment files
A list of paths to donor alignment (BAM/CRAM) files used in building the variant contexts.

### Donor variant files
A list of paths to variant (VCF/BCF) files used in building the variant contexts.

### Log file
Contains a log of steps and activities performed by VaSeBuilder. More info is added to the log when debug is activated.

### P-mode link file
When VaSeBuilder is run in P-mode, an extra file is created consisting of three columns

### Validation fastq files
Validation fastq files are valid fastq based on the provided template fastq files with acceptor and donor reads 
exchanged based on the variant contexts.

### Variant context file
The variant context file is a tab separated file containing the essential data of the variant contexts in 13 columns.
The context data is preceded by a header starting with a #.<br />

<u>Variant context file columns:</u>
* __ContextId:__ The identifier of the context (consists of the chromosome name and variant positions connected with 
an '_')
* __DonorSample:__ The sample name/identifier of the variant from which the context is constructed.
* __Chrom:__ The name of the chromosome the context is located on.
* __Origin:__ This is the genomic position of the donor variant the acceptor, donor and variant context 
* __Start:__ The starting, or leftmost genomic, position of the variant context.
* __End:__ The ending, or rightmost genomic, position of the variant context.
* __AcceptorContextLength:__ The length of the acceptor context that constructed the variant context.
* __DonorContextLength:__ The length of the donor context that constructed the variant context.
* __AcceptorReads:__ The number of acceptor reads and mates overlapping with the variant context.
* __DonorReads:__ The number of donor reads and mates overlapping with the variant context.
* __ADratio:__ The ratio of acceptor reads over donor reads that will be exchanged.
* __AcceptorReadIds:__ The read identifiers (one per mate pair) overlapping with the variant context.
* __DonorReadsIds:__ The read identifiers (one per mate pair) overlapping with the variant context.

The first line in each variant context file is the uuid of the VaSeBuilder run that constructed the variant contexts.
<br />

### Variant context statistics file

### Variant context bed file



## Debug output files
### Acceptor contexts file
### Acceptor context read positions file
### Acceptor context statistics file
### Acceptor context unmapped read mates file
### Donor contexts file
### Donor context read positions file
### Donor context statistics file
### Donor context unmapped read mates file
### Variant context acceptor read positions file
### Variant context donor read positions file
### Variant context unmapped acceptor read mates file
### Variant context unmapped donor read mate file
