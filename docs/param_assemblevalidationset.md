# AssembleValidationSet Parameters

##Required Inputs

- `-c / --varcon` __Variant context file:__ A variant context file containing the relevant information for the variant contexts represented by the spike-in BAMs and VCFs being used.
- `-kb / --spike-in-bam` __Donor spike-in BAM:__ Building block BAM files made using `BuildSpikeIns`. Multiple can be specified, space-separated.
	- `-kbL` Alternatively, a text file with a list of paths to BAM files, one per line.
- `-kv / --spike-in-vcf` __Donor spike-in VCF:__ Building block VCF files made using `BuildSpikeIns`. Multiple can be specified, space-separated.
	- `-kvL` Alternatively, a text file with a list of paths to VCF files, one per line.
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

* `-x / --hashtable` __Hashtable file__: Hashtable file of sample IDs of the samples provided in the run. Required if sample IDs are hashed in the provided variant context file. See [Hashing](nowhere).
* `-r / --reference` __Reference genome:__ A reference genome in FASTA format. The reference genome provided must be the same one used for alignment and variant calling of the acceptor sample and donor samples. If using CRAM input files, this parameter is required.
* `--seed` __Random seed:__ Integer used as the random seed to generate donor read positions. Donor reads are semi-randomly distributed across all template FastQ files to prevent alignment bias. This seed can be set to reproduce / shuffle validation sets.

---

##Output Controls

* `-o / --out-dir` __Output directory:__ Existing directory to write output files to. Default: Current directory
* `--fastq-out` __Output FastQ prefix:__ Prefix for output FastQ file names. Lane and pair number are appended automatically, such that a filename could be _MyPrefix\_L1\_R1.fastq_. Default: VaSe_<current_date\>
* `-l / --log` __Output log filename:__ Name for the log file that gets written during a VaSeBuilder run. Default: VaSeBuilder.log
* `--debug` __Debug logging:__ Output debug-level logging to log file and output additional information files. See [Debug output files](output_files.md#debug-output-files).

---

##Context Controls

* `--no-merge` __Disable same-sample context merging:__ When two variant contexts from the same sample overlap, the default behavior is to merge them into one wider context that includes both variants. Set to disable merging and default to any prioritization settings, if used.
* `--no-hash` __Disable sample ID hashing:__ Do not hash sample IDs when writing varcon files. See [Hashing](nowhere).

---
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


  --fastq-out <prefix> 
                        Prefix name for output FastQ files. Lane and pair number are appended automatically.
  --seed <int>          Random seed used to randomly distribute spike-in reads. (Default='VaSe_<date>'
  -av, --acceptor-vcf <vcf> 
                        Acceptor VCF file, used to make hybrid validation VCF.
  -c, --varcon <varconfile> 
                        Pre-made variant context file.
  -kb, --spike-in-bam <bam> [<bam2> ...] 
                        Pre-built spike-in BAM file(s).
  -kbL, --spike-in-bam-list <file> 
                        Pre-built spike-in BAM files listed per line in <file>.
  -kfq, --spike-in-fastq-list <file> 
                        Pre-built spike-in FastQ files with pairs listed tab-separated per line in <file>.
  -kv, --spike-in-vcf <vcf> [<vcf2> ...] 
                        Pre-built spike-in VCF file(s).
  -kvL, --spike-in-vcf-list <file> 
                        Pre-built spike-in VCF files listed per line in <file>.

Universal Options:
  -r, --reference <fasta> 
                        Reference sequence fasta
  -o, --out-dir <path> 
                        Output directory
  -l, --log <str>       Log output file name
  --debug               Log with maximum verbosity
