# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 13:29:17 2019

@author: medinatd
"""


import re
import sys
import pybedtools
# from dataclasses import dataclass

testfile = "/media/sf_ContinuousValidation/MVL/umcg_mvl_totaal_export_20190515.txt"
sort_order = [str(x) for x in range(1, 23)] + ["X", "Y", "MT"]


# =============================================================================
# @dataclass
# class MVL_record:
#     chromosome: str = "NULL"
#     start: str = "NULL"
#     stop: str = "NULL"
#     ref: str = "NULL"
#     alt: str = "NULL"
#     transcript: str = "NULL"
#     c_nomen: str = "NULL"
#     p_nomen: str = "NULL"
#     classification: str = "NULL"
#     variant_info: str = "NULL"
#     analysis_reference: str = "NULL"
#     comments: str = "NULL"
#     ID: str = None
#     ID_options: list = None
#     DNA: str = None
# =============================================================================


class MVL_record:
    def __init__(
            self, chromosome: str = "NULL", start: str = "NULL", stop: str = "NULL",
            ref: str = "NULL", alt: str = "NULL",
            transcript: str = "NULL", c_nomen: str = "NULL", p_nomen: str = "NULL",
            classification: str = "NULL", variant_info: str = "NULL",
            analysis_reference: str = "NULL", comments: str = "NULL",
            ID: str = None, ID_options: list = None, DNA: str = None,
            project: str = None, project_options: list = None
            ):
        self.chromosome = chromosome
        self.start = start
        self.stop = stop
        self.ref = ref
        self.alt = alt
        self.transcript = transcript
        self.c_nomen = c_nomen
        self.p_nomen = p_nomen
        self.classification = classification
        self.variant_info = variant_info
        self.analysis_reference = analysis_reference
        self.comments = comments
        self.ID = ID
        self.ID_options = ID_options
        self.DNA = DNA
        self.project = project
        self.project_options = project_options

    # XXX: Decide on DNA number comparison, short vs long.
    def __eq__(self, other):
        """Override the default Equals behavior."""
        if isinstance(other, self.__class__):
            return ([self.chromosome, self.start, self.stop, self.ref, self.alt, self.DNA] ==
                    [other.chromosome, other.start, other.stop, other.ref, other.alt, other.DNA])
        return False

    def __repr__(self):
        """Return string representation.

        Returns
        -------
        rep_string : str
        """
        rep_string = []
        for key in self.__dict__:
            rep_string.append(f"{key}=\'{self.__dict__[key]}\'")
        rep_string = "MVL_record:\n\t{}".format("\n\t".join(rep_string))
        return rep_string


class MVL:
    def __init__(self, MVL=None):
        self.header = []
        self.records = []
        self.chrom_list = []
        self.ref_fasta = None
        self.DNA_searched = (False, None)
        self.project_searched = False
        self.uniques_tagged = False
        if MVL is not None:
            self.read_from_file(MVL)
            self.set_chromosome_list()
        return

    def set_ref_fasta(self, reference):
        self.ref_fasta = reference
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

        reggies = [r"D(?:NA)?[0-9]{6}",
                   r"\w*?_D(?:NA)?[0-9]{6}_\w*",
                   r"(?<=A_)\w*?_DNA[0-9]{6}_\w*[0-9]",
                   r"(?<=A_)\w*?_DNA[0-9]{6}_[a-zA-Z0-9]*_[a-zA-Z0-9]*_[a-zA-Z0-9]*[0-9]"]
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

    def find_project_numbers(self):
        reg = re.compile(r"(QXT|QXTR|NGS|NGSR|CAR|5GPM)_?[0-9]{2,}", re.I)
        for record in self.records:
            regex_results = []
            regex_results.extend(re.findall(reg, record.analysis_reference))
            regex_results.extend(re.findall(reg, record.variant_info))
            regex_results.extend(re.findall(reg, record.comments))
            regex_results = list(set([res.replace("_", "") for res in regex_results]))

            num_results = len(regex_results)
            if num_results == 0:
                record.project = "NoneFound"
            elif num_results == 1:
                record.project = regex_results[0]
            elif num_results > 1:
                uniq_res = []
                for x in regex_results:
                    for y in regex_results:
                        if x in y and x != y:
                            break
                    else:
                        uniq_res.append(x)
                record.project = "Multiple"
                record.project_options = uniq_res

        self.project_searched = True

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
        self.uniques_tagged = True

    def subset_uniques(self):
        if not self.uniques_tagged:
            self.tag_uniques2()
        return [record for record in self.records
                if record.unique and record.ID not in ["NoneFound", "Multiple"]]

    def write_filter_list(self, outpath, short_ID=True, uniques_only=True):
        if uniques_only:
            if not self.uniques_tagged:
                self.tag_uniques2()
            choice = self.subset_uniques()
        else:
            choice = self.records
        with open(outpath, "w") as outfile:
            if short_ID:
                outfile.write("Sample\tChrom\tStart\tRef\tAlt\tClassification\tID\tProject\n")
                for record in choice:
                    outfile.write("\t".join([record.DNA, record.chromosome, record.start, record.ref, record.alt,
                                             record.classification, record.ID, record.project])
                                  + "\n")
            elif not short_ID:
                outfile.write("Sample\tChrom\tStart\tRef\tAlt\tClassification\n")
                for record in choice:
                    outfile.write("\t".join([record.ID, record.chromosome, record.start, record.ref, record.alt,
                                             record.classification])
                                  + "\n")

    def set_DNA_numbers(self):
        if not self.DNA_searched[0]:
            print("Run 'find_DNA_numbers' first.")
        for record in self.records:
            short_DNA = re.findall(r"D(?:NA)?[0-9]{6}", record.ID, re.IGNORECASE)
            if short_DNA:
                record.DNA = "DNA" + short_DNA[0][-6:]
            else:
                record.DNA = None

    @staticmethod
    def read_reference(reference, chrom, start, stop):
        """
        Return the nucleotides from positions in a reference fasta.

        Creates a temporary BED file and uses pybedtools (which wraps bedtools),
        to return the nucleotides at the BED positions from a provided reference
        fasta file.

        Parameters
        ----------
        reference : str
            Path to a reference fasta file.
        chrom : str
            Sequence (ie. chromosome) of the desired BED locus.
        start : int
            Start position of the desired locus.
        stop : int
            Stop position of the desired locus.

        Returns
        -------
        str
            Nucleotide sequence extracted from the reference.

        """
        bedpos = pybedtools.BedTool(f"{chrom}\t{start}\t{stop}\n",
                                    from_string=True)
        seq_out = bedpos.sequence(fi=reference)
        return seq_out.print_sequence().split("\n")[1]

    def realign_dot_notations(self):
        if self.ref_fasta is None:
            print("No reference provided.\nSet reference fasta first.")
        for record in self.records:
            if record.alt == ".":
                record.ref = MVL.read_reference(self.ref_fasta,
                                                record.chromosome,
                                                int(record.start) - 2,
                                                int(record.start) - 1) + record.ref
                record.alt = record.ref[0]
                record.start = str(int(record.start) - 1)
                record.stop = str(int(record.start) + len(record.ref) - 1)
            elif record.ref == ".":
                record.ref = MVL.read_reference(self.ref_fasta,
                                                record.chromosome,
                                                int(record.start) - 1,
                                                int(record.start))
                record.alt = record.ref + record.alt
                record.stop = int(record.start) + len(record.alt) - 1
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
    # myMVL = MVL(sys.argv[1])
    myMVL = MVL(testfile)
    myMVL.ref_fasta = "/media/sf_Ubuntu_Share/human_g1k_v37_phiX.fasta"
    for i in range(3, -1, -1):
        myMVL.find_DNA_numbers(i)
    myMVL.set_DNA_numbers()
    myMVL.tag_uniques2()
    # myMVL.find_DNA_numbers()
    # myMVL.add_loose_DNA_num()
