# Questions and Answers

__Q: Why is Python 3.7 or higher required?__  
A: When we started VaSeBuilder we wanted it to be ready for future Python versions. Furthermore, VaSeBuilder uses some Python code not available in earlier Python versions.

__Q: Why is pysam 0.15 or higher required?__  
A: 

__Q: What is the file command used for?__  
A: The linux file command is used to check that input files are of the required type. With alignment files for example, VaSeBuilder uses the file command to check that they are inded BAM or CRAM files rather than just examining the file extension.

__Q: Why does VASeBuilder have two specific steps?__  
A: We designed VaSeBuilder to have two specific steps to allow users to have some flexibility when creating .