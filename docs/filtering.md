# Filtering and Prioritization

When supplying an inclusion filter file, filtering and prioritization can be performed.

## Inclusion filter file
As specified in [inclusion filter file](input_file_formats.md#inclusion-filter-file), this user-created tab-separated text file specifies which variants to use from which sample while using `BuildSpikeIns` or `BuildValidationSet`. If an inclusion filter file is not provided, all variants from all VCFs provided will be used to create variant contexts. Otherwise, only VCF variants that match the included criteria will be used:

- __Sample:__ The VCF sample ID must match the specified sample
- __Chrom:__ The variant must be found on the specified chromosome
- __Pos:__ The variant's start position must match the specified position
- __Ref:__ The reference allele of the variant must match the specified reference sequence
- __Alt:__ The VCF variant record must contain the specified alternate allele, but can contain others as well.

These 5 criteria are the required columns and headers for the inclusion filter file:

```text
Sample	Chrom	Pos	Ref	Alt
sample1	3	12345678	A	T
sample2	4	1357911	GCC	G
sample3	7	24681012	T	TTT
```

Additional columns can be added, which can also be used as criteria for subsetting and prioritization:

```text
Sample	Chrom	Pos	Ref	Alt	Pathogenicity
sample1	3	12345678	A	T	benign
sample2	4	1357911	GCC	G	likely_pathogenic
sample3	7	24681012	T	TTT	benign
```


## Subsetting
Further criteria can be specified using `-s / --subset-filter`. Criteria are specified as `Column:value1,value2,...`, where `Column` refers to any column in the inclusion filter file, including custom columns. This setting further subsets by only allowing variants that additionally contain any of the specified `values` in `Column`, according to __OR__ logic. Multiple subset rules can be added, which then follow __AND__ logic, such that variants must satisfy each rule added for inclusion.

Example: `-s Pathogenicity:benign,likely_benign` will exclude the variant from sample2 in the above inclusion filter file, as it contains "likely\_pathogenic" in the Pathogenicity column, and not "benign" or "likely\_benign".

Example: `-s Pathogenicity:benign Chrom:4,5,6,7` will exclude the variants from sample1 and sample2. The variant from sample1 is listed as "benign" in the Pathogenicity column, but is not listed as being on chromosomes 4, 5, 6, or 7 in the Chrom column. The variant from sample2 is excluded again because it does not contain "benign" in the Pathogenicity column.

Example: `-s Flavor:chocolate` will not exclude any variants. Because there is no column called "Flavor", the criteria is ignored.

### Built-In Subset for Size
VaSeBuilder currently features a built-in subset for variant size that can be used without adding columns to the inclusion filter file. Using this subset follows the following format: `size:min[,max]`, and is based on the length of alleles reported in ref/alt.

Example inclusion filter file:

```text
Sample	Chrom	Pos	Ref	Alt	Pathogenicity
sample1	3	12345678	A	T	benign
sample2	4	1357911	GCC	G	likely_pathogenic
sample3	7	24681012	T	TTT	benign
sample4	X	11235813	CGCGCG	C	likely_benign
```

Example: `-s size:2` will only allow variants whose largest ref/alt allele includes &ge; 2 nucleotides. In the above inclusion filter file, the variant from sample1 will be excluded, as its longest allele has only 1 nucleotide.

Example: `-s size:0,3` will only include variants whose largest ref/alt allele is &ge; 0 nucleotides long, but &le; 3 nucleotides. In the above inclusion filter file, the variant from sample4 will be excluded, as it longest ref/alt allele is 6 nucleotides long.

## Prioritization


