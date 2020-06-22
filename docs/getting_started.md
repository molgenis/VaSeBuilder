# Getting started

## Get VaSeBuilder
VaSeBuilder can be obtained by downloading the [latest release](https://github.com/molgenis/VaSeBuilder/releases).
VaSeBuilder can also be easily installed via ```pip install vasebuilder```.

## Requirements
To run VaSeBuilder, the following software is required:  

* [Python 3.6 or higher](https://www.python.org)
* [pysam 0.15 or higher](https://pysam.readthedocs.io/en/latest/api.html)
* [numpy 0.18.1 or higher](https://numpy.org)
* [argon2-cffi 19.2.0 or higher](https://github.com/hynek/argon2-cffi)
* [pybedtools 0.8.1 or higher](http://daler.github.io/pybedtools/)
* [HTSLib 1.7 or higher](http://www.htslib.org)
* [SAMtools 1.7 or higher](http://www.htslib.org)
* Linux file command 5.32 or higher

## Run VaSeBuilder
To build a full validation set with VaSeBuilder run:
```
VaSeBuilder BuildValidationSet \  
  -b donor1.bam donor2.bam donor3.bam \
  -v donor1.vcf.gz donor2.vcf.gz donor3.vcf.gz \
  -a acceptor.bam \
  -1 acceptor_R1.fastq.gz -2 acceptor_R2.fastq.gz \
  -r reference_genome.fasta
```

or (when obtaining VaSeBuilder via the repository):
```
python -m VaSeBuilder
  -b donor1.bam donor2.bam donor3.bam \
  -v donor1.vcf.gz donor2.vcf.gz donor3.vcf.gz \
  -a acceptor.bam \
  -1 acceptor_R1.fastq.gz -2 acceptor_R2.fastq.gz \
  -r reference_genome.fasta
```

## Problem solving
