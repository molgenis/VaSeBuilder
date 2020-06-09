# Questions and Answers

__Q: Why is Python 3.7 or higher required?__  
A: When we started VaSeBuilder we wanted it to be ready for future Python versions. Furthermore, VaSeBuilder uses some Python code not available in earlier Python versions.

__Q: Why is pysam 0.15 or higher required?__  
A: VaSeBuilder uses pysam functionality that was added in version 0.l4.

__Q: What is the file command used for?__  
A: The linux file command is used to check that input files are of the required type. With alignment files for example, VaSeBuilder uses the file command to check that they are inded BAM or CRAM files rather than just examining the file extension.

__Q: Why does VASeBuilder have two specific steps?__  
A: We designed VaSeBuilder to have two specific steps to allow users to have some flexibility when creating .

__Q: Why are donor reads semi randomly added to the acceptor FastQ files?__  
A: Donor reads reads that need to be added are identified per variant context and therefore map near each other. If we would add these as is there would be 'blocks' in the FastQ files. This might introduce a potential bias during mapping as many consecutive reads would map to the same position when processing the validation set with a pipeline. Adding them randomly helps prevent this and also more closely resembles FastQ files from a sequencer.  
We chose to add them semi-randomly to ensure that if the exact same data, the same variant contexts and the same seed parameter is used, donor reads will be added at the same positions in the acceptor files.

__Q:__
