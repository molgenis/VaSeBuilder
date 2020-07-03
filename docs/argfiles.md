# Argument files

VaSeBuilder uses the Python Standard Library module __argparse__ for command line argument parsing. As such, it also supports argparse argument files. These are convenient, reusable sets of arguments that can be used instead of manually entering all arguments at the command line each time. They also provide a convenient way to track which commands were run with which arguments.

## Making an argument file
To write an argument file, add command line arguments to a text file, one per line, as they would appear on the command line beginning with the selected tool.

Example:

```text
BuildValidationSet  
--acceptor-bam GIAB_NA12878.bam  
--donor-vcf-list MyVCFList.txt  
--donor-bam-list MyBAMList.txt  
--reference hg37.fa  
--inclusion-filter MyFilterList.txt  
--acceptor-fq-r1 GIAB_L1_1.fq.gz GIAB_L2_1.fq.gz  
--acceptor-fq-r2 GIAB_L1_2.fq.gz GIAB_L2_2.fq.gz  
--fastq-out MyVaseSet  
--log MyLog.log  
#--no-hash
```

As shown in `#--no-hash`, individual lines can be commented out to conveniently disable/enable features between runs.

## Using an argument file

To use an argument file, run VaSeBuilder as follows, using `@` to denote the argument file:

```bash
VaSeBuilder @My_Argument_File.txt
```

In addition, individual arguments can be  added by specifying more arguments:

```bash
VaSeBuilder @My_Argument_File.txt --debug
```

Note that if an argument was already specified in the argument file, it will be overwritten if the argument is provided again on the command line. Arguments that accept multiple inputs, such as `--donor-bam`, will not have additional files appended, and will instead replace the entire argument.

Example:

```bash
VaSeBuilder @My_Argument_File.txt --inclusion-filter YourFilterList.txt
```

In the example above, the argument `--inclusion-filter YouFilterList.txt` will be used instead of the original argument `--inclusion-filter MyFilterList.txt` provided in the argument file.
