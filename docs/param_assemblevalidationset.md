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

* `-x / --hashtable` __Hashtable file__: Hashtable file of sample IDs of the samples provided in the run. See [Hashing](hashing.md).
* `-r / --reference` __Reference genome:__ A reference genome in FASTA format. The reference genome provided must be the same one used for alignment and variant calling of the acceptor sample and donor samples. If using CRAM input files, this parameter is required.
* `--seed` __Random seed:__ Integer used as the random seed to generate donor read positions. Donor reads are semi-randomly distributed across all template FastQ files to prevent alignment bias. This seed can be set to reproduce / shuffle validation sets.

---

##Output Controls

* `-o / --out-dir` __Output directory:__ Existing directory to write output files to. Default: Current directory
* `--fastq-out` __Output FastQ prefix:__ Prefix for output FastQ file names. Lane and pair number are appended automatically, such that a filename could be _MyPrefix\_L1\_R1.fastq_. Default: VaSe_<current_date\>
* `-l / --log` __Output log filename:__ Name for the log file that gets written during a VaSeBuilder run. Default: VaSeBuilder.log
* `--debug` __Debug logging:__ Output debug-level logging to log file and output additional information files. See [Debug output files](output_files.md#debug-output-files).
