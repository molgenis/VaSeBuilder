# Usage
VaSeBuilder has two overall steps: 1) building variant contexts, and 2) building validation FastQ files. These steps can be run as one continuous process, or as two separate processes. To run both steps in one go, runmode [```BuildValidationSet```](#vasebuilder-buildvalidationset) can be used. To run the steps separately, use [```BuildSpikeIns```](#buildspikeins), followed by [```AssembleValidationSet```](#assemblevalidationset).
VaSeBuilder can also be run with [```VaSeBuilder @argument_file.txt```](argfiles.md).

## BuildValidationSet
This run mode both builds variant contexts and builds a resulting set of FastQ files using these contexts. This run mode can, for example, be used when the donor data set is not large, when the first validation set is created, or when a small update to a validation set is required.  

_Example command:_

``` bash
VaSeBuilder BuildValidationSet \
    -b donor1.bam donor2.bam donor3.bam \
    -v donor1.vcf.gz donor2.vcf.gz donor3.vcf.gz \
    -f my_variants_of_interest.tsv \
    -a acceptor.bam \
    -1 acceptor_R1.fastq.gz -2 acceptor_R2.fastq.gz \
    -r reference_genome.fasta
```

---

## BuildSpikeIns
VaSeBuilder ```BuildSpikeIns``` establishes which donor and acceptor reads should be exchanged for each donor variant by establishing variant contexts. This first step therefore performs the essential work required to build the validation set. Users can run this step in different [output modes](#output-modes). This step always checks for overlapping variant context independent of the selected output mode.

To run this step, users need to at least supply donor data (alignment and variant files), an acceptor alignment file, and a genome reference. Other options such as filtering can also be set. Please see [Parameters](parameters.md) for more information.

_Example command:_

```bash
VaSeBuilder BuildSpikeIns \
    --output-mode P \
    -b donor1.bam donor2.bam donor3.bam \
    -v donor1.vcf.gz donor2.vcf.gz donor3.vcf.gz \
    -f my_variants_of_interest.tsv \
    -a acceptor.bam \
    -r reference_genome.fasta
```

### Output modes
```BuildSpikeIns``` can be run in different output modes. Each output mode outputs at least the [variant context file](output_files.md#variant-context-file) that contains all necessary information about the established variant contexts. Other output files are also written during this step, but which differ depending on selected output mode and whether VaSeBuilder is run in debug mode or not. Please see [Output files](output_files.md) for more information.

##### V-mode
V-mode only establishes variant contexts and outputs a variant context file. This mode differs from the other output modes in that it does not create alignment or variant output files. This mode can be helpful, for example, when you want to inspect established variant contexts using different input data, or to play with options without too producing too many output files.

##### P-mode
In addition to creating a variant context file, P-mode outputs an alignment and variant file for each used donor variant context. This mode thus creates small building blocks. Users can then select which buildings blocks ```AssembleValidationSet``` should use to build the validation set, allowing users to easily create different validation set 'flavours' using the same data.

##### D-mode
Like P-mode, D-mode outputs alignment and variant file buildings blocks, which users can later select to include in the validation set. Unlike P-mode, each alignment and variant file contains the combined data for all used donor variant contexts from an individual sample. This can be helpful when certain variants from a sample always need to be together. Due to the per-sample output files, the building blocks might be less flexible for creating different validation set flavours, but makes it easier to ensure that multiple variants from a single sample are included.  
_Note that D-mode is currently not yet available but will be implemented in the future._

##### A-mode
A-mode outputs only one alignment file and one variant file for all created variant contexts. 

---

## AssembleValidationSet
```AssembleValidationSet``` creates a validation set from a given set of building blocks. This step uses a set of acceptor FastQ files to create the validation set. The validation set itself therefore also consists of a set of FastQ files, equal to the number of acceptor FastQ files.  
Donor reads are semi randomly* added to the acceptor FastQ files, whilst acceptor reads within variant contexts are skipped.

_*When the exact same data and the same seed number is used, the donor reads should be inserted at the same position for reproducibility._

_Example command:_

```bash
VaSeBuilder AssembleValidationSet \
    -kb spike_in_1.bam spike_in_2.bam spike_in_3.bam \
    -kv spike_in_1.vcf.gz spike_in_2.vcf.gz spike_in_3.vcf.gz \
    -1 acceptor_R1.fastq.gz -2 acceptor_R2.fastq.gz \
    -c varcon_file.tsv
```
