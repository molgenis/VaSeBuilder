# Output files

## General output files

### Variant contexts

### Donors

### Others



## Debug output files
When VaSeBuilder is run in debug mode, additional output files are written to provide more information. These files can be helpful to get more information about the acceptor and donor contexts that help establish the variant contexts as well as more information about the variant context themselves.

### Acceptor and donor context files
First acceptor and donor context files are outputted. These files are very similar to the variant context file ```BuildSpikeIns``` produces and contain all essential data about the established acceptor and donor contexts. Since variant contexts are established by the acceptor and donor contexts, the left or right positions of these contexts might differ form that of the variant context for the same variant. Two files containing basic statistics about the acceptor and donor contexts are also written.

__Acceptor and donor context files written:__
* Acceptor contexts
* Donor contexts
* Acceptor context statistics
* Donor context statistics


### Unmapped read mate identifiers
When ```BuildSpikeIns``` establishes variant contexts, reads overlapping with the donor variant may have unmapped read mates. When VaSeBuilder is run in debug mode, the identifiers of reads that have an unmapped read mate are written to an output file for acceptor, donor and variant contexts. Note that there are two such files for variant contexts, one with acceptor reads and the other with donor reads. These files can have more read identifiers than the acceptor or donor contexts (Please see the 'Contexts' section for more information).

__Unmapped read mate files written:__
* Acceptor context unmapped read mate identifiers
* Donor context unmapped read mate identifiers
* Variant context unmapped acceptor read mate identifiers
* Variant context unmapped donor read mate identifiers


### Left and right read pair positions
Finally, files containing

__Left right position files written:__
* Acceptor context left and right most read pair positions
* Donor context left and right most read pair positions
* Variant context left and right most read pair positions
* Variant context left and right most raed pair positions
