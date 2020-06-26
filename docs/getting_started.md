# Getting started

## Get VaSeBuilder

VaSeBuilder can be obtained by downloading the [latest release](https://github.com/molgenis/VaSeBuilder/releases).
VaSeBuilder can also be easily installed via ```pip install vasebuilder```.

## Requirements
To run VaSeBuilder, the following software is required:  

* Libraries:
	* [HTSLib &ge; 1.7](http://www.htslib.org)
	* [SAMtools &ge; 1.7](http://www.htslib.org)
	* [Linux File(1) &ge; 5.32](https://github.com/file/file)
* Python:
	- [Python &ge; 3.6](https://www.python.org)
* Python packages:
	* [numpy &ge; 0.18.1](https://numpy.org)
	* [argon2-cffi &ge; 19.2.0](https://github.com/hynek/argon2-cffi)
	* [pysam &ge; 0.15](https://pysam.readthedocs.io/en/latest/api.html)
	* [pybedtools &ge; 0.8.1](http://daler.github.io/pybedtools/)

## Run VaSeBuilder
To build a full validation set with VaSeBuilder run:

```bash
VaSeBuilder BuildValidationSet \
  -b donor1.bam donor2.bam donor3.bam \
  -v donor1.vcf.gz donor2.vcf.gz donor3.vcf.gz \
  -a acceptor.bam \
  -1 acceptor_R1.fastq.gz -2 acceptor_R2.fastq.gz \
  -r reference_genome.fasta
```

or (when obtaining VaSeBuilder via the repository):

```bash
python -m VaSeBuilder \
  -b donor1.bam donor2.bam donor3.bam \
  -v donor1.vcf.gz donor2.vcf.gz donor3.vcf.gz \
  -a acceptor.bam \
  -1 acceptor_R1.fastq.gz -2 acceptor_R2.fastq.gz \
  -r reference_genome.fasta
```

## Problem solving
