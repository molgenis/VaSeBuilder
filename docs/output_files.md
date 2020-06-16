# Output files

## General output files

### Variant contexts
For each run of ```BuildSpikeIns``` or ```BuildValidationSet``` VaSeBuilder writes multiple output files. First three output files with information about variant contexts are written. The variant contexts file contains information about the created variant contexts, such as the window size and reads in the context. the variant context statistics file contains basic statistics about the variant contexts. The third file contains the variant contexts in BED format, which makes it easy to check the variant context locations and windows in a genome browser.

__Variant context files written:__

* [Variant contexts](output_file_formats.md#variant-context-file)
* [Variant context statistics](output_file_formats.md#variant-context-statistics-file)
* [Variant contexts BED](output_file_formats.md#variant-context-bed-file)

### Donor lists
The second set of files are two files that list the used donor alignment and variant files respectively, to construct the variant contexts.

__Donor list files written:__

* [Donor alignment file list](output_file_formats.md#donor-alignment-files)
* [Donor variant file list](output_file_formats.md#donor-variant-files)

### FastQ files
VaSeBuilder also output a set of FastQ files when ```AssembleValidationSet``` or ```BuildValidationSet``` is run. These FastQ files will have acceptor reads in variant contexts exchanged by donor reads in variant contexts. The number of FastQ files is equal to the number of acceptor FastQ files supplied.

__FastQ files written:__

* R1 FastQ file(s)
* R2 FastQ files(s)

### Others
Currently VaSeBuilder also outputs two other files: the log file and the hashtable. Note that the hashtable is not written if sample ID hashing has been turned off.


## Debug output files
When VaSeBuilder is run in debug mode, additional output files are written to provide more information. These files can be helpful to get more information about the acceptor and donor contexts that help establish the variant contexts as well as more information about the variant context themselves.

### Acceptor and donor context files
First acceptor and donor context files are outputted. These files are very similar to the variant context file ```BuildSpikeIns``` produces and contain all essential data about the established acceptor and donor contexts. Since variant contexts are established by the acceptor and donor contexts, the left or right positions of these contexts might differ form that of the variant context for the same variant. Two files containing basic statistics about the acceptor and donor contexts are also written.

__Acceptor and donor context files written:__

* [Acceptor contexts](output_file_formats.md#acceptordonor-contexts-file)
* [Donor contexts](output_file_formats.md#acceptordonor-contexts-file)
* [Acceptor context statistics](output_file_formats.md#acceptordonor-context-statistics-file)
* [Donor context statistics](output_file_formats.md#acceptordonor-context-statistics-file)


### Unmapped read mate identifiers
When [```BuildSpikeIns```](usage.md#buildspikeins) establishes [variant contexts](contexts.md#variant-context), reads overlapping with the donor variant may have unmapped read mates. When VaSeBuilder is run in debug mode, the identifiers of reads that have an unmapped read mate are written to an output file for [acceptor](output_file_formats.md#acceptordonor-context-unmapped-read-mates-file), [donor](output_file_formats.md#acceptordonor-context-unmapped-read-mates-file) and [variant contexts](output_file_formats.md#variant-context-unmapped-acceptordonor-read-mates-file). Note that there are two such files for variant contexts, one with acceptor reads and the other with donor reads. These files can have more read identifiers than the acceptor or donor contexts (Please see the 'Contexts' section for more information).

__Unmapped read mate files written:__

* [Acceptor context unmapped read mate identifiers](output_file_formats.md#acceptordonor-context-unmapped-read-mates-file)
* [Donor context unmapped read mate identifiers](output_file_formats.md#acceptordonor-context-unmapped-read-mates-file)
* [Variant context unmapped acceptor read mate identifiers](output_file_formats.md#variant-context-unmapped-acceptordonor-read-mates-file)
* [Variant context unmapped donor read mate identifiers](output_file_formats.md#variant-context-unmapped-acceptordonor-read-mates-file)


### Left and right read pair positions
Finally, files containing the left and right most positions for each read pair per variant context can be written. Users can use these files to make start/end position plots that can show outliers.

__Left right position files written:__

* Acceptor context left and right most read pair positions
* Donor context left and right most read pair positions
* Variant context left and right most read pair positions
* Variant context left and right most read pair positions
