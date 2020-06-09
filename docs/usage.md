# Usage
VaSeBuilder offers users to run the program in different ways, known as runmodes. VaSeBuilder can be run via ```python vase.py```, followed by a runmode and setting the required general parameters as well as other required and optional parameters. There are three runmodes: ```BuildValidationSet```, ```BuildSpikeIns``` and ```AssembleValidationSet```.

## VaSeBuilder steps
VaSeBuilder has two specific steps that can be run together or separately. To run both steps in one go, VaSeBuilder ```BuildValidationSet``` can be used. Users can however also run both steps separately using VaSeBuilder ```BuildSpikeIns``` and VaSeBuilder ```AssembleValidationSet```.

### BuildSpikeIns
VaSeBuilder ```BuildSpikeIns``` is the step that establishes which donor and acceptor reads should be exchanged for each donor variant. This first step therefore performs the essential work required to build the validation set. Users can run this step in different output modes (further described below).
To run this step users need at least supply donor data (alignment and variant files), an acceptor alignment file and a genome reference. Other options such as filtering can also be set (Please see section 'Parameters' for more information)


#### Output modes
VaSeBuilder ```BuildSpikeIns``` can be run in different output modes. Each output mode outputs at least the variant context file that contains all necessary information about the established variant contexts (see 'Output files' for more information). Other output files are also written during this step, but which differs per selected output mode and whether VaSeBuilder is run in debug mode or not.

##### V-mode
V-mode only establishes variant contexts and outputs a variant context file. This mode therefore differs from the other output modes, like P-mode in that it does not create alignment and variant output files. This mode can be helpful for example when you want to inspect established variant contexts using different input data or options without too much output files.

##### P-mode
P-mode outputs an alignment and variant file for each used donor variant. It therefore does not check for overlapping variant contexts as each variant is processed individually. This mode thus creates small building blocks. Users can then select which buildings blocks VaSeBuilder ```AssembleValidationSet``` should use to build the validation set, allowing users to easily create different valdiation set 'flavours' using the same data.

_Example command:_
```
python vase.py BuildSpikeIns \
    --output-mode P \
    -b donor1.bam donor2.bam donor3.bam \
    -v donor1.vcf.gz donor2.vcf.gz donor3.vcf.gz \
    -f my_variants_of_interest.tsv \
    -a acceptor.bam \
    -r reference_genome.fasta
```

##### D-mode
D-mode outputs alignment and variant files per sample. This mode creates per sample buildings blocks users can later select to include in the validation set. This can be helpful when certain variants from a sample always need to be together. Due to the per sample output files, the building blocks might be less flexible for in creating different validation set flavours, it makes it easier to include multiple variants from a single sample.  
_Note that D-mode is currently not available yet but will bne implemented in the future._

##### A-mode
A-mode outputs one alignment and one variant file for all created variant contexts. 


### AssembleValidationSet
VaSeBuilder ```AssembleValidationSet``` is the step that creates the validation set.
