# Hashing

VaSeBuilder was designed for use with clinical patient data, making privacy and security a concern. To increase anonymization of data, VaSeBuilder by default uses the Argon2 hashing algorithm to hash the sample ID's found in VCF and BAM files.

This has implications for `BuildSpikeIns` and `BuildValidationSet`. Both tools produce a variant context file, and each each variant context has its original sample ID recorded. If hashing is turned on (default), these sample ID's are replaced with a hashed value.

Sample IDs are also replaced in the building blocks produced by `BuildSpikeIns`:

* BAM files contain sample IDs in the `@RG` header line `SM` field. These are replaced with the hashed ID. Currently, the `LB` field is also replaced with the sample ID hash, though the `ID` field is not altered.
* VCF files contain sample IDs in the samples' respective column headers.
	* When using `--runmode P`, each VCF sample column header is replaced with its sample ID hash.
	* When using `--runmode A`, the single VCF output file has only one sample column, whose ID is replaced with the placeholder "VaSeBuilder".
	* Additionally, some annotation tools may output sample IDs in `##INFO` header line descriptions. Currently, VCF files output by VaSeBuilder remove all `##INFO` description fields to avoid sample ID leaks.

## Hash Table
Argon2 is a one-way hashing algorithm, meaning there is no way to directly retrieve the original sample ID from the hash. However, it may be necessary to re-access the original files from which validation sets were made, particularly when reusing variant context files.

To make this possible, VaSeBuilder produces a hash table, which is simply a table listing the sample IDs and their hashes. This obviously circumvents anonymization, and as such this file should only be accessible wherever it is safe to keep the private data itself.

Sample ID hashes are random and are not reproducible, such that the same sample will receive different hashed IDs if run through VaSeBuilder multiple times. However, if a hash for the sample ID exists in a VaSeBuilder hash table, the hash table can be provided to reuse this hash, instead of generating a new one.

If hashing was enabled while using `BuildSpikeIns` with `--runmode V`, the produced hash table must be supplied to reuse the variant context file so that VaSeBuilder can connect the hashed sample IDs to their original files.