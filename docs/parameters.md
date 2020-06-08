# VaSeBuilder Parameters

## Global parameters
These parameters are always required indepent of the selected output mode.

### General parameters
* __[-m / --output-mode] Selected output mode:_ This option allows users to select which output mode VaSeBuilder should be run in. The output mode can be specified with a single letter with A (A-mode), D (D-mode), P (P-mode) and V (V-mode) as accepted values. Note that D-mode has not yet been implemented.
* __[-r / --reference] Genome reference:__ One single reference can be provided and should be in FASTA format. Furthermore, this genome reference needs to be the reference used to process (read mapping, variant calling, etc) both the acceptor sample and donor sammples.
* __[-o / --out-dir] Output directory:__ Path to an existing directory where VaSebuilder should write the output files to.
* __[-l / --log] Log file:__ Users can provide a name for the log file that gets written during a VaSebuilder run. The default name for a log file is 'VaSeBuilder.log'. When multiple VaSeBuilder runs are performed this option can be useful to differentiate log files.

### Alignment file parameters
There are two mutually exclusive parameters for providing alignment (BAM/CRAM) files. Note that alignment files can also consist of a mixture of BAM and CRAM files.
* __[-b / --donor-bam] Donor alignment files:__ Donor sample alignment files can be provided with each filepath separated by a ','. Each sample should have only one file provided. Provided donor files that do not exist are skipped in subsequent analysis steps.
* __[-bL / --donor-bam-list] Donor alignment list file:__ Donor sample alignment files are provided by means of a list file. This list file needs to have each single alignment file on a separate line. Donor files that do not exist are skipped in subsequent analysis steps.

### Variant file parameters
As in the case with the donor alignment files above, there are also two mutually exclusive parameters for providing variant (VCF/bgzipped VCF) files.
* __[-v / --donor-vcf] Donor variant files:__ Donor sample variant files can be provided with each filepath separated by a ','. Each sample should have only one file provided. Provided donor files that do not exist are skipped.
* __[-vL / --donor-variant-list] Donor variant list file:__ Donor sample variant files are provided by means of a list file. As with the alignment list file, each variant file should be on a separate line. Donor files that do not exist are skipped.

It is important to note that for each sample both an alignment file and variant file are required in order to be used by VaSeBuilder for context creation. Samples with an alignment or variant file will therefore be excluded.


### Acceptor FastQ file parameters
* __[-1 / --acceptor-fq-r1] Acceptor FastQ R1:__ R1 acceptor FastQ files to use as acceptor/template files, to spike donor reads into, can be provided with each filepath separated by a ','.
* __[-1L / --acceptor-fq-r1-list] Acceptor FastQ R1 list:__ R1 acceptor FastQ file can also be provided by means of a list file. Each line should have only one FastQ file.
* __[-2 / --acceptor-fq-r2] Acceptor FastQ R2:__ R2 acceptor FastQ files to use as acceptor/template files, to spike donor reads into, can be provided with each filepath separated by a ','.
* __[-2L / --acceptor-fq-r2-list] Acceptor FastQ R2 list:__ R2 acceptor FastQ file can also be provided by means of a list file. Each line should have only one FastQ file.
* __[--fastq-out] FastQ out name:__ When VaSeBuilder outputs a set of FastQ files with spiked in variants, the default output name 'VaSe_' followed by the data and either R1 or R2 and lane number like L1 or L2, etc. Users can can specify a prefix that will replace 'VaSe_' and the date. The R1/R2 and lane numbers are added after the prefix.
* __[--seed] Random seed:__ Integer to set the seed to semi-randomly distributed donor reads over the template FastQ files. This is in order to prevent donor reads that map to same location from forming blocks in the FastQ file.
* __[-av / --acceptor-vcf] Acceptor VCF:__ 


## Optional parameters

### Filtering parameters
* __[-f / --inclusion-filter] Inclusion filter file:__ Users can provide a filter file specifying which variants to use (include). The filter file requires at least several required columns after which optional custom columns can be added. Please see input formats for more information. This parameter needs to be set in order to be able to use the subset-filter and prioritization filter arguments.
* __[-s / --subset-filter] Subset filter:__ 
* __[-p / --prioritization] Priority filter:__ 

### Context creation
* __[--no-merge] Don't merge:__ Overlapping variant contexts from the same sample are not merged and saved as two separate variants contexts. This option is useful when running VaSeBuilder in P-mode as both variant contexts will be saved and written.
* __[-vo / -- varcon-out] Variant context file outname:__ By default, the variant context output file is named 'VaSe_date' with date being the current date. Users can specifiy a name for the output file.

### Miscellaneous
* __[--no-hash] No sample ID hashing:__ By default VaSebuilder uses Argeon2 to hash sample identifiers. This flag allows users to disable this behaviour for the current run.
* __[-x / --hashtable] Alreay existing hashtable:__ Users can also provide an already created (Argon2) hashtable to hash sample identifiers with.


## Runmode specific parameters
### A-mode
### D-mode
### P-mode
### V-mode


## Future parameters
VaSeBuilder also has a few parameters that have not yet been implemented but we aim to do in the future.

* __[-m D / --output-mode D] VaSeBuilder running in D-Mode:__ 
* __[--suppress-conflit-check] Suppress checking for conflicts between variant contexts:__ Conflicts caused between two variant contexts from different samples that overlap are ignored. This option can be helpful when VaSeBuilder is run in P-mode.
* __[--add-secondary-variants] Add secondary variants to VCF output file:__ If a variant excluded by filtering overlaps with an included variant, the excluded variant will still be written to the included variant VCF output file. This option can be helpful when the excluded variant is from the same sample as the included variant as the excluded variant is on the same set of reads. Variant callers will therefore most likely still call the excluded variant. This allows 'suprise' variants to show up.
* __[-O / --output-type] Type of output file to write:__ Users can specify which type of output the validation set with exchanged reads and therefore spikedin variants should be. These output types will be BAM ('B'), UBAM ('U') or a set of FASTQ files ('F'). If not set, the default output type will be a BAM file.
