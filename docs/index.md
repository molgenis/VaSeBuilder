# VaSeBuilder
## About VaSeBuilder
Often when a computational analysis consisting of multiple steps is run routinely, a pipeline is constructed. Pipelines connect the individual analysis steps and automate their execution and the transitions to the next steps. Pipelines should be evaluated to determine whether each step and transition works correctly. Simultaneously, validation of the pipeline to determine whether it can indeed provide relevant results is equally important. Pipelines can be validated in separate and isolated runs, but this requires a lot of time. With VaSeBuilder we aim to perform pipeline validation and sample analyses simultaneously.

### General workings
VaSeBuilder (Validation Set Builder) can be used to construct a set of FastQ files containing reads with selected and therefore known variants. To construct the set of FastQ files, alignment and variant data from samples (the donors) are used to modify a template (the acceptor), which could be NA12878 for example.
VaSeBuilder first collects donor read pairs with the desired variants from provided samples. Template read pairs overlapping with the variant position are then exchanged for sample reads. Donor read pairs are collected from samples by identifying aligned read pairs (from BAM/CRAM files) of which at least one read overlaps with the variant position. Acceptor reads are collected in the same way from the acceptor alignment (BAM/CRAM) file.

To determine which reads should be exchanged, multiple windows called [Contexts](contexts.md) are constructed.

### Important to know
For the best results, please use sample data and an acceptor that have been sequenced using the same preparation kit, capturing kit, and sequencing platform.

VaSeBuilder works best with whole-exome sequencing (WES) data.


## To get started
To get started with VaSeBuilder run <command to place here>
