# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 13:29:17 2019

@author: medinatd
"""


import re
import sys
from dataclasses import dataclass

testfile = "../../MVL/umcg_mvl_totaal_export_20190515.txt"
sort_order = [str(x) for x in range(1, 23)] + ["X", "Y", "MT"]


@dataclass
class MVL_record:
    chromosome: str = "NULL"
    start: str = "NULL"
    stop: str = "NULL"
    ref: str = "NULL"
    alt: str = "NULL"
    transcript: str = "NULL"
    c_nomen: str = "NULL"
    p_nomen: str = "NULL"
    classification: str = "NULL"
    variant_info: str = "NULL"
    analysis_reference: str = "NULL"
    comments: str = "NULL"
    ID: str = None
    ID_options: list = None
    DNA: str = None

    # XXX: Decide on DNA number comparison, short vs long.
    def __eq__(self, other):
        """Override the default Equals behavior"""
        if isinstance(other, self.__class__):
            return ([self.chromosome, self.start, self.stop, self.ref, self.alt, self.DNA] ==
                    [other.chromosome, other.start, other.stop, other.ref, other.alt, other.DNA])
        return False


class MVL:
    def __init__(self, MVL=None):
        self.header = []
        self.records = []
        self.chrom_list = []
        self.DNA_searched = (False, None)
        self.uniques_tagged = False
        if MVL is not None:
            self.read_from_file(MVL)
            self.set_chromosome_list()
        return

    def read_from_file(self, MVL):
        with open(MVL, errors="replace") as MVL_in:
            header = MVL_in.readline()
            rawMVL = MVL_in.readlines()
        self.header = header.strip().split("\t")
        rawMVL = [line.strip().split("\t") for line in rawMVL]
        self.records = [MVL_record(**dict(zip(self.header, record))) for record in rawMVL]
        return

    def set_chromosome_list(self):
        self.chrom_list = set([x.chromosome for x in self.records])
        return

    def find_DNA_numbers(self, strictness=2, replace=False):
        # if sys.version_info[1] == 8:
        #     self.__findpy38()
        #     return

        reggies = [r"\w*D(?:NA)?[0-9]{6}\w*",
                   r"(?<=A_)\w*?_D(?:NA)?[0-9]{6}_\w*",
                   r"\w*?_D(?:NA)?[0-9]{6}_\w*"]
        reg = re.compile(reggies[strictness], re.I)

        for record in self.records:
            if not replace and record.ID is not None and record.ID != "NoneFound":
                continue

            regex_results = []
            regex_results.extend(re.findall(reg, record.analysis_reference))
            regex_results.extend(re.findall(reg, record.variant_info))
            regex_results.extend(re.findall(reg, record.comments))
            regex_results = list(set([res.lstrip("A_") for res in regex_results]))

            num_results = len(regex_results)
            if num_results == 0:
                record.ID = "NoneFound"
            elif num_results == 1:
                record.ID = regex_results[0]
            elif num_results > 1:
                uniq_res = []
                for x in regex_results:
                    for y in regex_results:
                        if x in y and x != y:
                            break
                    else:
                        uniq_res.append(x)
                record.ID = "Multiple"
                record.ID_options = uniq_res
        self.DNA_searched = (True, strictness)

    def sort_by_genomic_position(self):
        extra_chroms = [x for x in self.chrom_list if x not in sort_order]
        extra_chroms.sort()
        sort_order.extend(extra_chroms)
        self.records.sort(key=lambda x: int(x.start))
        self.records.sort(key=lambda x: sort_order.index(x.chromosome))

    def tag_uniques(self, designate_primaries=False):
        for chrom in self.chrom_list:
            chrom_subset = [x for x in self.records if x.chromosome == chrom]
            for record in chrom_subset:
                if chrom_subset.count(record) == 1:
                    record.unique = True
                else:
                    record.unique = False
        self.uniques_tagged = True

    def tag_uniques2(self):
        for chrom in self.chrom_list:
            chrom_unique_list = []
            chrom_subset = [x for x in self.records if x.chromosome == chrom]
            for record in chrom_subset:
                if chrom_unique_list.count(record) == 0:
                    chrom_unique_list.append(record)
                    record.unique = True
                else:
                    record.unique = False

    def subset_uniques(self):
        if not self.uniques_tagged:
            self.tag_uniques()
        return [record for record in self.records
                if record.unique and record.ID not in ["NoneFound", "Multiple"]]

    def write_filter_list(self, outpath, uniques_only=True):
        if uniques_only:
            if not self.uniques_tagged:
                self.tag_uniques()
            choice = self.subset_uniques()
        else:
            choice = self.records
        with open(outpath, "w") as outfile:
            outfile.write(f"Sample\tChrom\tPos\tRef\tAlt\n")
            for record in choice:
                outfile.write("\t".join([record.ID, record.chromosome, record.start, record.ref, record.alt]) + "\n")

    def add_loose_DNA_num(self):
        if not self.DNA_searched[0]:
            print("Run 'find_DNA_numbers' first.")
        for record in self.records:
            short_DNA = re.findall(r"D(?:NA)?[0-9]{6}", record.ID, re.IGNORECASE)
            if short_DNA:
                record.DNA = "DNA" + short_DNA[0][-6:]
            else:
                record.DNA = None

# =============================================================================
#             if (re.findall(reg, record.analysis_reference)):
#                 r = re.findall(reg, record.analysis_reference)
#                 self._assign_ID(record, r)
#
#             elif (re.findall(reg, record.variant_info)):
#                 r = re.findall(reg, record.variant_info)
#                 if len(set(r)) > 1:
#                     record.ID = "Multiple"
#                 else:
#                     record.ID = r[0]
#
#             elif (re.findall(reg, record.comments)):
#                 r = re.findall(reg, record.comments)
#                 if len(set(r)) > 1:
#                     record.ID = "Multiple"
#                 else:
#                     record.ID = r[0]
#
#             else:
#                 record.ID = "NULL"
#             record.ID = record.ID.strip("A_")
# =============================================================================

# =============================================================================
#     def __findpy38(self):
#         reg = re.compile("(?<=A_)\w*?_D(?:NA)?[0-9]{6}_\w*", re.I)
#         for record in self.records:
#             if (r := re.findall(reg, record.analysis_reference)):
#                 record.ID = r[0]
#             elif (r := re.findall(reg, record.variant_info)):
#                 if len(set(r)) > 1:
#                     record.ID = "Multiple"
#                 else:
#                     record.ID = r[0]
#             elif (r := re.findall(reg, record.comments)):
#                 if len(set(r)) > 1:
#                     record.ID = "Multiple"
#                 else:
#                     record.ID = r[0]
#             else:
#                 record.ID = "NULL"
#             record.ID.strip("A_")
# =============================================================================


if __name__ == "__main__":
    myMVL = MVL(sys.argv[1])
    # myMVL = MVL(testfile)
    # myMVL.find_DNA_numbers()
    # myMVL.add_loose_DNA_num()
