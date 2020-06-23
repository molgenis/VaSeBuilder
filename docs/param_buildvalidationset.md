# BuildValidationSet Parameters

##Required Inputs

- `-b / --donor-bam` __Donor alignment files:__ Indexed BAM or CRAM files for donor samples. Multiple can be specified, space-separated.
	- `-bL` Alternatively, a text file with a list of paths to BAM or CRAM files, one per line.
- `-v / --donor-vcf` __Donor variant files:__ Indexed VCF.GZ or BCF files for donor samples. Multiple can be specified, space-separated.
	- `-vL` Alternatively, a text file with a list of paths to VCF.GZ or BCF files, one per line.
- `-a / --acceptor-bam` __Acceptor alignment file:__ Acceptor BAM or CRAM file to use as the template when making variant contexts.
- `-1 / --acceptor-fq-r1` __Acceptor FastQ R1 files:__ Multiple can be specified, space-separated.
	- `-1L` Alternatively, a text file with a list of paths to FastQ R1 files, one per line.
- `-2 / --acceptor-fq-r2` __Acceptor FastQ R2 files:__ Multiple can be specified, space-separated.
	- `-2L` Alternatively, a text file with a list of paths to FastQ R2 files, one per line.

_Note that acceptor FastQ R1 and R2 files must be provided in the order of their pairing. Example:_

```bash
-1 acceptor_L1_R1.fq.gz acceptor_L2_R1.fq.gz \
-2 acceptor_L1_R2.fq.gz acceptor_L2_R2.fq.gz
```

---

##Optional Inputs

* `-x / --hashtable` __Hashtable file__: If a VaSeBuilder hashtable has already been created for samples in this run, their hashes can be reused by specifying the file. See [Hashing](nowhere).
* `-r / --reference` __Reference genome:__ A reference genome in FASTA format. The reference genome provided must be the same one used for alignment and variant calling of the acceptor sample and donor samples. If using CRAM input files, this parameter is required.
* `--seed` __Random seed:__ Integer used as the random seed to generate donor read positions. Donor reads are semi-randomly distributed across all template FastQ files to prevent alignment bias. This seed can be set to reproduce / shuffle validation sets.

---

##Output Controls

* `-o / --out-dir` __Output directory:__ Existing directory to write output files to. Default: Current directory
* `--fastq-out` __Output FastQ prefix:__ Prefix for output FastQ file names. Lane and pair number are appended automatically, such that a filename could be _MyPrefix\_L1\_R1.fastq_. Default: VaSe_<current_date\>
* `-vo / --varcon-out` __Output variant context filename:__ Filename for the variant context output file. Default: VaSe_<current_date\>.varcon
* `-l / --log` __Output log filename:__ Name for the log file that gets written during a VaSeBuilder run. Default: VaSeBuilder.log
* `--debug` __Debug logging:__ Output debug-level logging to log file and output additional information files. See [Debug output files](output_files.md#debug-output-files).

---

##Context Controls

* `--no-merge` __Disable same-sample context merging:__ When two variant contexts from the same sample overlap, the default behavior is to merge them into one wider context that includes both variants. Set to disable merging and default to any prioritization settings, if used.
* `--no-hash` __Disable sample ID hashing:__ Do not hash sample IDs when writing varcon files. See [Hashing](nowhere).

---

##Filtering Controls

* `-f / --inclusion-filter` __Inclusion filter:__ A tab-separated text file with a list of variants to include, formatted with the following mandatory columns with column headers:
	- _Sample_: Sample ID in which the variant occurs. Sample ID must match the sample ID in its corresponding BAM and VCF files.
	- _Chrom_: Chromosome (or contig) in which the variant occurs, as indicated in its VCF record.
	- _Pos_: Genomic position at which the variant starts, as indicated in its VCF record.
	- _Ref_: Reference allele for the variant, as indicated in its VCF record.
	- _Alt_: Alternative allele for the variant, as indicated in its VCF record. Does not need to match ALL alternative alleles, as this is mostly used for concordance checking.

Additional custom columns may be added with unique column headers for use in further filtering using the `-s` and `-p` flags.

* `-s / --subset-filter` __Inclusion filter subset:__ Subset the provided inclusion filter file by specified criteria. Subsets are specified as `Column:value1,value2,...`, where `Column` refers to a column in the inclusion filter file, and `value` refers to a desired value found in that column. Any variant containing any of the specified values in the specified column will be used, while others will be ignored. User-specified custom columns beyond the 5 required columns can be used as well. Multiple subsetting rules can be specified, space-separated.
	* Example: Specifying `-s Chrom:1,3,X Ref:A,T` will keep only variants recorded as being on chromosomes 1, 3, or X with a reference allele recorded as an A or T, while other variants are discarded.
* `-p / --prioritization` __Prioritization rules:__ In the event that two variant contexts overlap and cannot be merged (see [Merging](nowhere)), only one context is kept while the other is discarded. When considering which context to discard, prioritization rules can be specified to preferentially keep higher-priority contexts based on the properties of their included variants as recorded in the inclusion filter. Priorities are specified as `Column:value1,value2,...`, where `Column` refers to a column in the inclusion filter file, and `value` refers to a desired value found in that column. Priorities are set from highest to lowest, such that a variant with `value1` in `Column` is preferentially kept over any variant with `value2`, and both are preferentially kept over any variant with any other value in `Column`. User-specified custom columns beyond the 5 required columns can be used. Multiple prioritization rules can be specified, space-separated, and are also considered highest to lowest.
	* Example, Specifying `--prioritization Pathogenicity:P,LP,VOUS Type:InDel`:
		* Context 1 contains a variant with "P" and "SNP" recorded in the user-specified "Pathogenicity" and "Type" columns of an inclusion filter file
		* Context2 contains a variant with "B" and "InDel" recorded for the same columns
		* Result: Context 1 is kept and context 2 is discarded, as "Pathogenicity" takes priority over "Type", and within "Pathogenicity", "P" takes priority over "B"

See [Filtering and Prioritization](nowhere) for more information.

---

## BuildSpikeIns

## AssembleValidationSet

### General parameters
* `-m / --output-mode` __Selected output mode:__ This option allows users to select which output mode VaSeBuilder should be run in. The output mode can be specified with a single letter with A (A-mode), D (D-mode), P (P-mode) and V (V-mode) as accepted values. Note that D-mode has not yet been implemented.
* `-r / --reference` __Genome reference:__ One single reference can be provided and should be in FASTA format. Furthermore, this genome reference needs to be the reference used to process (read mapping, variant calling, etc) both the acceptor sample and donor samples.
* `-o / --out-dir` __Output directory:__ Path to an existing directory where VaSeBuilder should write the output files to.
* `-l / --log` __Log file:__ Name for the log file that gets written during a VaSeBuilder run. The default name for a log file is 'VaSeBuilder.log'. When multiple VaSeBuilder runs are performed this option can be useful to differentiate log files.

### Alignment file parameters
There are two mutually exclusive parameters for providing alignment (BAM/CRAM) files. Note that alignment files can also consist of a mixture of BAM and CRAM files.

* `-b / --donor-bam` __Donor alignment files:__ Donor sample alignment files can be provided with each file path separated by a space. Users should supply only one alignment file per sample.
* `-bL / --donor-bam-list` __Donor alignment list file:__ Donor sample alignment files are provided by means of a list file. This list file needs to have each single alignment file on a separate line.

### Variant file parameters
As in the case with the donor alignment files above, there are also two mutually exclusive parameters for providing variant (VCF/bgzipped VCF) files.

* `-v / --donor-vcf` __Donor variant files:__ Donor sample variant files can be provided with each file path separated by a space. Users should supply only one variant file per sample.
* `-vL / --donor-variant-list` __Donor variant list file:__ Donor sample variant files are provided by means of a list file. As with the alignment list file, each variant file should be on a separate line.

It is important to note that for each sample both an alignment file and variant file are required in order to be used by VaSeBuilder for context creation. Samples with an alignment or variant file will therefore be excluded.


### Acceptor FastQ file parameters
* `-1 / --acceptor-fq-r1` __Acceptor FastQ R1:__ R1 acceptor FastQ files to use as acceptor/template files, to spike donor reads into, can be provided with each file path separated by a space.
* `-1L / --acceptor-fq-r1-list` __Acceptor FastQ R1 list:__ R1 acceptor FastQ file can also be provided by means of a list file. Each line should have only one FastQ file.
* `-2 / --acceptor-fq-r2` __Acceptor FastQ R2:__ R2 acceptor FastQ files to use as acceptor/template files, to spike donor reads into, can be provided with each file path separated by a space.
* `-2L / --acceptor-fq-r2-list` __Acceptor FastQ R2 list:__ R2 acceptor FastQ file can also be provided by means of a list file. Each line should have only one FastQ file.
* `--fastq-out` __FastQ out name:__ When VaSeBuilder outputs a set of FastQ files with spiked in variants, the default output name 'VaSe_' followed by the data and either R1 or R2 and lane number like L1 or L2, etc. Users can can specify a prefix that will replace 'VaSe_' and the date. The R1/R2 and lane numbers are added after the prefix.
* `--seed` __Random seed:__ Integer to set the seed to semi-randomly distributed donor reads over the template FastQ files. This is in order to prevent donor reads that map to same location from forming blocks in the FastQ file.
* `-av / --acceptor-vcf` __Acceptor VCF:__ 


## Optional parameters

### Filtering parameters
* `-f / --inclusion-filter` __Inclusion filter file:__ Path to a filter file specifying which variants to use (include). The filter file requires at least several required columns after which optional custom columns can be added. Please see input formats for more information. This parameter needs to be set in order to be able to use the ```--subset-filter``` and ```--prioritization``` arguments.
* `-s / --subset-filter` __Subset filter:__ Column names and values as criteria for the inclusion filter. Entries not satisfying one of the values will be excluded. A column and values to use as filter criteria can be specified as ```--subet-filter <column name>: <value> <value>```.
* `-p / --prioritization` __Variant priority:__ Column names and values to use as prioritization. Prioritization assigns priorities to column from left to right with the first column obtaining the highest priority, the second slightly lower, etc. The same applies to the provided values for each column. Values in the prioritization column(s) that were not mentioned will be assigned the lowest priority. Prioritization can be specified as ```--prioritization <column name>: <value> <value>```

### Context creation
* `--no-merge` __Don't merge:__ Overlapping variant contexts from the same sample are not merged and saved as two separate variants contexts. This option is useful when running VaSeBuilder in P-mode as both variant contexts will be saved and written.
* `-vo / -- varcon-out` __Variant context file out name:__ By default, the variant context output file is named 'VaSe_date' with date being the current date. Users can specify a name for the output file.

### Miscellaneous
* `--no-hash` __No sample ID hashing:__ By default VaSeBuilder uses Argon2 to hash sample identifiers. This flag allows users to disable this behaviour for the current run.
* `-x / --hashtable` __Already existing hashtable:__ Path to an already created (Argon2) hashtable to hash sample identifiers with.



## Future parameters
VaSeBuilder also has a few parameters that have not yet been implemented but we aim to do in the future.

* `-m D / --output-mode D` __VaSeBuilder running in D-Mode:__ 
* `--suppress-conflict-check` __Suppress checking for conflicts between variant contexts:__ Conflicts caused between two variant contexts from different samples that overlap are ignored. This option can be helpful when VaSeBuilder is run in P-mode.
* `--add-secondary-variants` __Add secondary variants to VCF output file:__ If a variant excluded by filtering overlaps with an included variant, the excluded variant will still be written to the included variant VCF output file. This option can be helpful when the excluded variant is from the same sample as the included variant as the excluded variant is on the same set of reads. Variant callers will therefore most likely still call the excluded variant. This allows 'surprise' variants to show up.
* `-O / --output-type` __Type of output file to write:__ Type of output the validation set with exchanged reads and therefore spikedin variants should be. These output types will be BAM ('B'), UBAM ('U') or a set of FastQ files ('F'). If not set, the default output type will be a BAM file.
