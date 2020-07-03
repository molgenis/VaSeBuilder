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

* `-x / --hashtable` __Hashtable file__: If a VaSeBuilder hashtable has already been created for samples in this run, their hashes can be reused by specifying the file. See [Hashing](hashing.md).
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
* `--no-hash` __Disable sample ID hashing:__ Do not hash sample IDs when writing varcon files. See [Hashing](hashing.md).

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
* `-p / --prioritization` __Prioritization rules:__ In the event that two variant contexts overlap and cannot be merged (see [Merging](contexts.md#merging)), only one context is kept while the other is discarded. When considering which context to discard, prioritization rules can be specified to preferentially keep higher-priority contexts based on the properties of their included variants as recorded in the inclusion filter. Priorities are specified as `Column:value1,value2,...`, where `Column` refers to a column in the inclusion filter file, and `value` refers to a desired value found in that column. Priorities are set from highest to lowest, such that a variant with `value1` in `Column` is preferentially kept over any variant with `value2`, and both are preferentially kept over any variant with any other value in `Column`. User-specified custom columns beyond the 5 required columns can be used. Multiple prioritization rules can be specified, space-separated, and are also considered highest to lowest.
	* Example, Specifying `--prioritization Pathogenicity:P,LP,VOUS Type:InDel`:
		* Context 1 contains a variant with "P" and "SNP" recorded in the user-specified "Pathogenicity" and "Type" columns of an inclusion filter file
		* Context2 contains a variant with "B" and "InDel" recorded for the same columns
		* Result: Context 1 is kept and context 2 is discarded, as "Pathogenicity" takes priority over "Type", and within "Pathogenicity", "P" takes priority over "B"

See [Filtering and Prioritization](nowhere) for more information.
