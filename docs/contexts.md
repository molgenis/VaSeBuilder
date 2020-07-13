# Contexts
During the process of exchanging acceptor reads with donor reads, VaSeBuilder creates multiple windows containing reads we call Contexts. For each variant, overlapping read pairs are collected. The size of the context window for each variant is determined by the left and rightmost genomic position of the overlapping read pairs. Sometimes mates from reads overlapping with the variant can be mapped far away from the variant location. Such reads are considered outliers and are filtered by their start/stop position using Tukey's Fences method when constructing the window sizes. These outlier reads are however still spiked into the acceptor.
<br /><br />
![Context creation](img/vsb.gif)

## Types of contexts
### Acceptor & donor contexts
Acceptor and donor contexts are the windows constructed by acceptor and donor read pairs, respectively, overlapping with the variant position. The windows of an acceptor and donor context for a variant might differ in terms of length, leftmost and rightmost genomic positions due to read mapping differences. Important to note is that only read pairs with at least one read overlapping with the variant are used to construct the acceptor and donor context.

### Variant context
A Variant context for each variant is constructed by determining the leftmost and rightmost genomic positions from both the acceptor and donor context. Once this window has been established, acceptor and donor reads within this window are collected and exchanged (if the variant context is part of a Super context, than all acceptor and donor reads overlapping with the super context are exchanged.)

### Super context
Super contexts consist of two or more variant contexts that have been merged into a single context, and therefore contain multiple variants. Super contexts can occur when a sample has multiple variants that have overlapping variant contexts.

---

## Merging
When multiple variant contexts overlap, one of two things can occur:

- One context is discarded, while the other is kept
- The contexts are merged into one larger super context

This behavior depends on a several factors:

- __Sample of origin:__ Variant contexts from different samples cannot be merged.
- __User-defined behavior:__ While the default behavior is to merge overlapping variant contexts from the same sample, users can specify `--no-merge` to disable this.
- __User-defined prioritizations:__ Users can specify variant prioritization behaviors by using `--inclusion-filter` and `--prioritization`. The overlapping variant context containing the highest-priority variants is retained, while the other is discarded.
- __Order:__ If merging is not possible / disabled, and if there is a tie in priority or no prioritizations are set, then whichever variant context was made first during the run is retained, while the newer context is discarded.
