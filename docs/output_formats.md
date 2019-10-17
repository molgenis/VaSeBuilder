# VaSeBuilder output file formats

## Standard output files
VaSeBuilder runs output several output files containing the essential data about the run. Some files are only outputted 
when running specific modes.<br /><br />



### Donor alignment files
A list of paths to donor alignment (BAM/CRAM) files used in building the variant contexts. A list could be:
<br /><br />
_/path/to/donor1.bam<br />
/path/to/donor2.cram<br />
/path/to/donor3.bam_
<br /><br />


### Donor insert positions
VaSeBuilder add donor reads at semi random positions in the resulting validation set. The recorded position is the 
position of a full read (in a fastq a donor read with insert position 1 would be written to lines 5,6,7 and 8).
<br />

_\#VBUUID {vbuuid}<br />
ReadId &emsp;R1_InsertPos &emsp;FastqR1Out &emsp;R2_InsertPos &emsp;FastqR2Out<br />
dReadId1 &emsp;15 &emsp;/testdata/vase_R1.fastq &emsp; 15 &emsp;/testdata/vase_R2.fastq_
<br /><br />


### Donor variant files
A list of paths to donor variant (VCF/BCF) files used in building the variant contexts. A list could be:<br /><br />
_/path/to/donor1.vcf<br />
/path/to/donor2.bcf<br />
/path/to/donor3.vcf_
<br /><br />


### Log file
Contains a log of steps and activities performed by VaSeBuilder. More info is added to the log when debug is activated.
The log file contains four fields separated by a tab: Date and time, name of the logger (VaSe_Logger), log level, 
message. The first two lines of the log always display the issued command to run VaSeBuilder and the VaSeBuilder uuid 
and creation date and time A VaSeBuilder log file therefore follows the format:<br />

_YYYY-MM-DD &emsp;HH:MM:SS,SSS &emsp;INFO &emsp;python vase.py -c test.cfg<br />
YYYY-MM-DD &emsp;HH:MM:SS,SSS &emsp;INFO &emsp;VaSeBuilder: {uuid} ; YYYY-MM-DD HH:MM:SS,SSS<br />
YYYY-MM-DD &emsp;HH:MM:SS,SSS &emsp;INFO &emsp;Running VaSeBuilder in F-mode_
<br /><br />


### P-mode link file
When VaSeBuilder is run in P-mode, donor fastq files are outputted for each variant context. To link the variant 
context file and the outputted donor fastq files, an extra file is created consisting of three columns:
<br /><br />
_contextid_1 &emsp;/path/to/fastq_R1.fq &emsp;/path/to/fastq_R1.fq<br />
contextid_2 &emsp;/path/to/fastq_R2.fq &emsp;/path/to/fastq_R2.fq_
<br />

Each line represents a variant context. The first line in the P-mode link file is the uuid of the VaSeBuilder run. This 
is used to the variant context file and the fastq files for a certain P-mode run.
<br /><br />


### Validation fastq files
Validation fastq files are valid fastq files based on the provided template fastq files with acceptor and donor reads 
exchanged based on the variant contexts. Currently (September 24th, 2019) donor reads are added at the end of the 
validation fastq files.
<br /><br />


### Variant context file
The variant context file is a tab separated file containing the essential data of the variant contexts in 13 columns. 
This file is one of the important output files as it contains the windows created to search reads and which acceptor 
reads are exchanged for which donor reads. The context data is preceded by a header starting with a #. The default 
output name for the variant context file is _varcon.txt_.<br />

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

The first line in each variant context file is the uuid of the VaSeBuilder run that constructed the variant contexts. A 
variant context file looks like:<br />

_#VBUUID: {uuid}<br />
\#ContextId &emsp;DonorSample &emsp;Chrom &emsp;Origin &emsp;Start &emsp;End &emsp;AcceptorContextLength 
&emsp;DonorContextLength &emsp;AcceptorReads &emsp;DonorReads &emsp;ADratio &emsp;AcceptorReadsIds 
&emsp;DonorReadIds<br />
1_100 &emsp;Sample1 &emsp;1 &emsp;100 &emsp;50 &emsp;150 &emsp;75 &emsp;75 &emsp;2 &emsp;2 &emsp;1.0 
&emsp;aReadId1,aReadId2 &emsp;dReadId1,dReadId2_
<br /><br />


### Variant context statistics file
The variant context statistics file contains some statistics about the reads and mates overlapping with the variant 
contexts. Per variant context the read statistics consists, in order, of the average and median read length, q-score 
and mapq values for acceptor and donor reads overlapping with the variant context.
<br /><br />
_#ContextId &emsp;Avg_ALen &emsp;Avg_DLen &emsp;Med_ALen &emsp;Med_DLen &emsp;Avg_AQual &emsp;Avg_DQual &emsp;Med_AQual 
&emsp;Med_DQual &emsp;Avg_AMapQ &emsp;Avg_DMapQ &emsp;Med_AMapQ &emsp;Med_DMapQ<br />
1_100 &emsp;151 &emsp;151 &emsp;151.0 &emsp;151.0 &emsp;36.5 &emsp;35.6 &emsp;37.0 &emsp;38.0 &emsp;54.6 &emsp;55.7 
&emsp;60.0 &emsp;60.0_
<br /><br />


### Variant context bed file
A BED file with the variant context entries in four columns: Chromosome name, context start, context end, context 
identifier. The output name for the BED file is _variantcontexts.bed_. The resulting file looks like:<br />

_1 &emsp;100 &emsp;1000 &emsp;1_500<br />
2 &emsp;200 &emsp;2000 &emsp;2_1000_
<br /><br />


## Debug output files
When debug is set to True multiple extra files are outputted as well providing more in depth information.
<br /><br />


### Variant context acceptor/donor read positions file
Much like the acceptor/donor read positions the variant context acceptor positions and variant context donor positions 
files save the leftmost R1 and rightmost R2 acceptor and donor read positions respectively.
<br /><br />


### Variant context unmapped acceptor/donor read mates file
Records the identifiers of reads that have an unmapped read mate in a tab separated manner. The file contains three 
columns: Context identifier, sample name/identifier, read identifiers.<br /><br />

_#ContextId &emsp;SampleId &emsp;ReadIds<br />
1_100&emsp;Sample1 &emsp;uReadId1,uReadId2,uReadId3_


### Acceptor/Donor contexts file
The acceptor/donor context file is essentially a variant context file but for acceptor or donor contexts respectively. 
The acceptor/donor context file differs in the number of columns and only has two additional columns after ContextId - 
End, namely NumOfReads and ReadIds.

_VBUUID: {uuid}<br />
\#ContextId &emsp;DonorSample &emsp;Chrom &emsp;Origin &emsp;Start &emsp;End &emsp;AcceptorContextLength 
&emsp;DonorContextLength &emsp;NumOfReads &emsp;ReadIds<br />
1_100 &emsp;Sample1 &emsp;1 &emsp;100 &emsp;50 &emsp;150 &emsp;2 &emsp;dReadId1,dReadId2_
<br /><br />


### Acceptor/Donor context statistics file
Similar to the variant context statistics file, the acceptor/donor context statistics file contains basic statistics of 
the reads overlapping with acceptor or donor contexts.
<br /><br />
_#ContextId &empsp;Avg_ReadLen &emsp;MedReadLen &emsp;AvgReadQual &emsp;Med_ReadQual &emsp;AvgReadMapQ 
&emsp;MedReadMapQ<br />
1_100 &emsp;151 &emsp;151.0 &emsp;37.9 &emsp;39.0 &emsp;57.1 &emsp;60.0_
<br /><br />


### Acceptor/Donor context read positions file
The read positions contains the leftmost genomic positions of all R1 reads and rightmost genomic positions of all R2 
reads in all acceptor or donor contexts. These positions can be plotted, if needed, to gain more insight into the outer 
ranges of the reads associated with a variant context. The file contains three columns: context identifier, left 
positions, right positions. Each leftmost and rightmost position is separated by a comma.<br /><br />
_#ContextId &emsp;LeftPos &emsp;RightPos<br />
1_100 &emsp;25,50,75 &emsp;75,50,25_
<br /><br />


### Acceptor/Donor context unmapped read mates file
Reads overlapping with the acceptor/donor contexts with unmapped read mates encountered during a VaSeBuilder run are 
written the the acceptor_unmapped.txt or donor_unmapped.txt
<br /><br />
