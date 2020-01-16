# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:51:25 2020

@author: tdmedina
"""

import pysam


def combine_validation_variants(vcf_path_file, outpath, sample="Prometheus"):
    with open(vcf_path_file) as list_file:
        vcf_paths = list_file.readlines()
    vcf_paths = [x.strip() for x in vcf_paths]
    variant_list = []
    for vcf_path in vcf_paths:
        vcf = pysam.VariantFile(vcf_path)
        variant_list.extend(var for var in vcf.fetch())
        vcf.close()

    variant_list.sort(key=lambda x: int(x.pos))
    sort_order = [str(x) for x in range(1, 23)] + ["X", "Y", "MT"]
    variant_list.sort(key=lambda x: sort_order.index(x.chrom))

    fields = ["fileformat", "filter", "alt", "format", "contig", "reference", "info"]
    # header_list = [var.header.records for var in variant_list]
    header_list = [str(x) for var in variant_list
                   for x in var.header.records
                   if str(x).lstrip("#").split("=")[0].lower() in fields]
    print(header_list)
    header_list_unique = []
    for header in header_list:
        if header in header_list_unique:
            continue
        header_list_unique.append(header)
    header_list_unique.append(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
    with open(outpath, "w") as outfile:
        outfile.writelines(header_list_unique)
        for var in variant_list:
            outfile.write(var.__str__())


if __name__ == "__main__":
    import sys
    combine_validation_variants(sys.argv[1], sys.argv[2])
