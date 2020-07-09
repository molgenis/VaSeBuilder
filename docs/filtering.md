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

* Example: `-s Pathogenicity:benign,likely_benign` will exclude the variant from sample2 in the above inclusion filter file, as it contains "likely\_pathogenic" in the Pathogenicity column, and not "benign" or "likely\_benign".

* Example: `-s Pathogenicity:benign Chrom:4,5,6,7` will exclude the variants from sample1 and sample2. The variant from sample1 is listed as "benign" in the Pathogenicity column, but is not listed as being on chromosomes 4, 5, 6, or 7 in the Chrom column. The variant from sample2 is excluded again because it does not contain "benign" in the Pathogenicity column.

* Example: `-s Flavor:chocolate` will not exclude any variants. Because there is no column called "Flavor", the criteria is ignored.

### Built-In Subset for Size
VaSeBuilder currently features a built-in subset for variant size that can be used without adding columns to the inclusion filter file. Using this subset follows the format `size:min[,max]`, and is based on the length of the largest allele reported in ref/alt.

Example inclusion filter file:

```text
Sample	Chrom	Pos	Ref	Alt	Pathogenicity
sample1	3	12345678	A	T	benign
sample2	4	1357911	GCC	G	likely_pathogenic
sample3	7	24681012	T	TTT	benign
sample4	X	11235813	CGCGCG	C	likely_benign
```

* Example: `-s size:2` will only allow variants whose largest ref/alt allele includes &ge; 2 nucleotides. In the above inclusion filter file, the variant from sample1 will be excluded, as its longest allele has only 1 nucleotide.

* Example: `-s size:0,3` will only include variants whose largest ref/alt allele is &ge; 0 nucleotides long, but &le; 3 nucleotides. In the above inclusion filter file, the variant from sample4 will be excluded, as it longest ref/alt allele is 6 nucleotides long.

## Prioritization
In the event that two variant contexts overlap and cannot be merged (see [Merging](contexts.md#merging)), only one context is kept while the other is discarded. When considering which context to discard, prioritization rules can be specified to preferentially keep higher-priority contexts based on the properties of their included variants as recorded in the inclusion filter. Priorities are specified as `Column:value1,value2,...`, where `Column` refers to a column in the inclusion filter file, and `value` refers to a desired value found in that column. Priorities are set from highest to lowest, such that a variant with `value1` in `Column` is preferentially kept over any variant with `value2`, and both are preferentially kept over any variant with any other value in `Column`. User-specified custom columns beyond the 5 required columns can be used. Multiple prioritization rules can be specified, space-separated, and are also considered highest to lowest.

Example inclusion filter file:

```text
Sample	Chrom	Pos	Ref	Alt	Pathogenicity	Effect
sample1	3	12345678	A	T	benign synonymous
sample2	3	12345678	A	AG	pathogenic	frameshift_insertion
sample3	3	12345678	A	AGGG	VOUS	inframe_insertion
sample4	3	12345678	A	ATAG	pathogenic	premature_stop
```

The four variants above occur at the same locus in three different samples. There will be an overlap in their variant contexts, since they are located at the same locus, but they will not be merged because they come from different samples.

* Example: No prioritization set. Whichever variant is processed first will be kept, while the other three variants will be discarded when their overlaps are detected.

* Example: Specifying `--prioritization Pathogenicity:pathogenic,likely_pathogenic,VOUS`. Because `benign` is not listed in the prioritization, the variant from sample1 will be assigned the lowest priority and will be discarded, since the variants from the other samples will take priority. Because `pathogenic` is listed before `VOUS` in the prioritization settings, the context from sample3 will be discarded, as the variants from sample2 and sample4 will take precedence. Finally, because no other criteria are set, the decision between the variants from samples 2 and 4 will be made based on whichever variant was processed first.

* Example: Specifying `--prioritization Pathogenicity:pathogenic Effect:premature_stop `. In this case, `Effect` is listed after `Pathogenicity`, so comparisons based on `Effect` values will only be used to break ties. Comparisons will proceed as in the above example; however, because the final two variants from sample2 and sample4 are tied based on `Pathogenicity`, their `Effect` values will be compared, and the variant from sample4 will take precedence.

### Built-In Prioritization for Size
VaSeBuilder currently features a built-in prioritization for variant size that can be used without adding columns to the inclusion filter file. To specify this prioritization, use either `size:greater` or `size:less`. Contexts will then be prioritized based on the lengths of their longest alleles, with `greater` prioritizing longer-allele variants, and `less` prioritizing smaller-allele variants.

### Super Context Prioritization Behavior
In the event that a [super context](contexts.md#super-context) (the result of merging two overlapping contexts from the same sample) is created, its prioritization values are created by taking the highest individual values from each of its component contexts.

Example inclusion filter file:

```text
Sample	Chrom	Pos	Ref	Alt	Pathogenicity	TiTv
sample1	5	24681012	A	G	benign transition
sample1	5	24681015	T	A	pathogenic	transversion
sample2	5	24681012	A	T	pathogenic	transversion
```

* Example: `--prioritization Pathogenicity:pathogenic TiTv:transition` In the above example, the two variant contexts from sample1 will be merged (if enabled) into one super context. The variant context created for the variant from sample2, however, cannot be merged because it is from a separate sample. When comparing contexts to determine which should be kept, the combined super context will now have `pathogenic` and `transition` as its values, which are the highest priority values from each of its component variants. As such, the combined super context will be kept over the context created for sample2, since `transition` takes priority over `transversion`.