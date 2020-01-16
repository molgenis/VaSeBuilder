# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 10:26:25 2019

@author: medinatd
"""

# _default_annotation_order = ["Classification", "Ref", "Alt"]


class Context:
    def __init__(self, ContextId, DonorSample, Chrom, Origin, Start, End,
                 AcceptorContextLength, DonorContextLength,
                 AcceptorReads, DonorReads, ADratio,
                 AcceptorReadsIds, DonorReadIds, group=None,
                 **kwargs):
        self.ContextId = ContextId
        self.DonorSample = DonorSample
        self.Chrom = Chrom
        self.Origin = Origin
        self.Start = Start
        self.End = End
        self.AcceptorContextLength = AcceptorContextLength
        self.DonorContextLength = DonorContextLength
        self.AcceptorReads = AcceptorReads
        self.DonorReads = DonorReads
        self.ADratio = ADratio
        self.AcceptorReadsIds = AcceptorReadsIds
        self.DonorReadIds = DonorReadIds
        self.group = group
        self.Annotations = kwargs

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.ContextId, self.DonorSample) == (other.ContextId, other.DonorSample)
        return False

    def to_str(self):
        out_str = "\t".join([self.ContextId, self.DonorSample, self.Chrom,
                             self.Origin, self.Start, self.End,
                             self.AcceptorContextLength, self.DonorContextLength,
                             self.AcceptorReads, self.DonorReads, self.ADratio,
                             self.AcceptorReadsIds, self.DonorReadIds, str(self.group),
                             "\t".join(self.Annotations.values())])
        return out_str

class Combined_Varcon:
    def __init__(self):
        self.contexts = []
        self.conflicts = []
        self.grouped = False

    def add_contexts_from_file(self, infile):
        with open(infile, "r") as varcon_in:
            new_varcons = varcon_in.readlines()
        header = new_varcons[1].strip().strip("#").split()
        new_varcons = [varcon.strip().split("\t") for varcon in new_varcons if not varcon.startswith("#")]
        for varcon in new_varcons:
            self.contexts.append(Context(**dict(zip(header, varcon))))

    def add_contexts_from_list(self, infile):
        with open(infile, "r") as list_in:
            paths = list_in.readlines()
        paths = [path.strip() for path in paths]
        for path in paths:
            self.add_contexts_from_file(path)

    def build_checklist(self, context_list):
        checklist = {str(chrom): [] for chrom in list(range(1, 23)) + ["X", "Y"]}
        for record in context_list:
            checklist[record.Chrom].append(record)
        return checklist

    def divide_conflicts(self):
        new_contexts = []
        new_conflicts = []

        checklist = self.build_checklist(self.contexts)

        for variant in self.contexts:
            for check in checklist[variant.Chrom]:
                # if (variant.ContextId == check.ContextId
                #         and variant.DonorSample == check.DonorSample):
                #     continue
                if variant == check:
                    continue
                elif (int(variant.Start) <= int(check.End)
                        and int(check.Start) <= int(variant.End)):
                    new_conflicts.append(variant)
                    break
            else:
                new_contexts.append(variant)

        self.contexts = new_contexts
        self.conflicts = new_conflicts

    def display_conflicts(self):
        for x in self.conflicts:
            print(x)

    def display_counts(self):
        print(f"Contexts: {len(self.contexts)}")
        print(f"Conflicts: {len(self.conflicts)}")

    def write_contexts(self, outfile="combined_contexts.varcon"):
        writer = ["\t".join([str(y) for y in [x.__dict__.values()][:-1]])
                  + "\t".join(x.Annotations.values())
                  + "\n"
                  for x in self.contexts]
        with open(outfile, "w") as contexts_out:
            contexts_out.write("#VBUUID: Conflict_Free\n")
            contexts_out.write("#ContextId\tDonorSample\tChrom\tOrigin\tStart\tEnd\t"
                               "AcceptorContextLength\tDonorContextLength\t"
                               "AcceptorReads\tDonorReads\tADratio\t"
                               "AcceptorReadsIds\tDonorReadIds\t"
                               + "\t".join(self.contexts[0].Annotations.keys())
                               + "\n")
            contexts_out.writelines(writer)

    def write_contexts2(self, outfile, contexts):
        writer = [context.to_str() + "\n" for context in contexts]
        header = ("#ContextId\tDonorSample\tChrom\tOrigin\tStart\tEnd\t"
                  "AcceptorContextLength\tDonorContextLength\t"
                  "AcceptorReads\tDonorReads\tADratio\t"
                  "AcceptorReadsIds\tDonorReadIds\tConflictGroup\t"
                  + "\t".join(contexts[0].Annotations.keys())
                  + "\n")
        with open(outfile, "w") as contexts_out:
            contexts_out.write(f"#VBUUID: {outfile}\n{header}")
            contexts_out.writelines(writer)

    def write_conflicts(self, outfile="conflicts.varcon"):
        writer = ["\t".join([str(y) for y in [x.__dict__.values()][:-1]])
                  + "\t".join(x.Annotations.values())
                  + "\n"
                  for x in self.conflicts]
        with open(outfile, "w") as contexts_out:
            contexts_out.write("#VBUUID: Conflict_Only\n")
            contexts_out.write("#ContextId\tDonorSample\tChrom\tOrigin\tStart\tEnd\t"
                               "AcceptorContextLength\tDonorContextLength\t"
                               "AcceptorReads\tDonorReads\tADratio\t"
                               "AcceptorReadsIds\tDonorReadIds\t"
                               + "\t".join(self.conflicts[0].Annotations.keys())
                               + "\n")
            contexts_out.writelines(writer)
# =============================================================================
#     def write_varcon(self, list_in, outfile):
# 
#         if list_in[0].Annotations is not None:
#             for i in range(len(list_in[0].Annotations))
#             header = header + "\tAnnotation"
# =============================================================================
    def group_overlaps(self):
        i = 0
        for x in self.conflicts:
            x.group = []
        for x in self.conflicts:
            added = False
            for y in [context for context in self.conflicts if context.Chrom == x.Chrom]:
                if x == y:
                    continue
                elif (int(x.Start) <= int(y.End)
                        and int(y.Start) <= int(x.End)):
                    if (set(x.group) & set(y.group)):
                        continue
                    else:
                        x.group.append(i)
                        y.group.append(i)
                        added = True
                        continue
            if added:
                i += 1
        self.grouped = True

        collapsed_multis = self._collapse_complex_groups()
        for x in self.conflicts:
            for collapsed_multi in collapsed_multis:
                if set(x.group).issubset(collapsed_multi):
                    x.group = min(collapsed_multi)
                    break
            else:
                x.group = x.group[0]
        return

# =============================================================================
#         complex_groups = [x.group for x in self.conflicts if len(x.group) > 1]
#         compound_complex = []
#         for complex1 in complex_groups:
#             new_complex = set(complex1)
#             for complex2 in complex_groups:
#                 if set(new_complex) & set(complex2):
#                     new_complex = set(complex1) | set(complex2)
#             for complex_2 in complex_groups:
#                 if set(complex1) & set(complex2):
#                     new_complex = set(complex1) | set(complex2)
#             re_base_group = []
#             for conflict in self.conflicts:
#                 if conflict.
# =============================================================================
    def _collapse_complex_groups(self):
        multi_groups = [x.group for x in self.conflicts if len(x.group) > 1]
        multi_groups.reverse()
        combined_multis = []
        while multi_groups:
            combined_multi = set(multi_groups.pop())
            for other_multi in multi_groups:
                if combined_multi & set(other_multi):
                    combined_multi = combined_multi | set(other_multi)
            for other_multi in list(combined_multis):
                if combined_multi & other_multi:
                    combined_multi = combined_multi | other_multi
                    del combined_multis[combined_multis.index(other_multi)]
            combined_multis.append(combined_multi)
        return combined_multis

    def interactive_resolver(self):
        group_list = sorted(list(set([x.group for x in self.conflicts])))
        # group_list.sort(lambda x: "," in x)
        # group_list_complex = set([int(x) for y in group_list for x in y.split(",") if "," in y])
        # group_list_simple = [x for x in group_list
        #                      if "," not in group_list
        #                      and x not in ]
        self.keepers = []
        for i in group_list:
            conflict_group = [x for x in self.conflicts if x.group == i]
            num = 1
            print(f"\nConflict group: {i} / {group_list[-1]}")
            header = "Num\tContextID\tAConLen\tDConLen\tAReads\tDReads\tADratio"
            for head in self.conflicts[0].Annotations.keys():
                header += f"\t{head:20}"
            # print("Num\tAConLen\tDConLen\tAReads\tDReads\tADratio" + "\t".join(self.conflicts[0].Annotations.keys())
            print(header)
            for x in conflict_group:
                info = (f"{num}:\t{x.ContextId}\t{x.AcceptorContextLength}\t{x.DonorContextLength}\t"
                        f"{x.AcceptorReads}\t{x.DonorReads}\t{x.ADratio:.6}")
                for annot in x.Annotations.values():
                    if len(annot) > 25:
                        info += f"\t{annot:.25}...({len(annot)} total)"
                    else:
                        info += f"\t{annot:20}"
                # print(f"{num}:\t{x.ContextId}\t{x.AcceptorContextLength}\t{x.DonorContextLength}\t"
                #       f"{x.AcceptorReads}\t{x.DonorReads}\t{x.ADratio:.6}\t"
                #       + "\t".join(x.Annotations.values()))
                print(info)
                num += 1
            chosen = False
            while not chosen:
                choices = input()
                try:
                    choices = list(set([int(x) for x in input().split(",")]))
                    for choice in choices:
                        self.keepers.append(conflict_group[choice - 1])
                    chosen = True
                except (IndexError, ValueError):
                    print("Invalid choice. Try again.")

            # conflict_group[choice - 1].group = f"{conflict_group[choice - 1].group}_resolved"
            self.keepers.append(conflict_group[choice - 1])
        print("End of simple conflicts.")

# =============================================================================
#     def prioritize_pathogenicity(self):
#         self.resolved = []
#         group = []
#         for context in self.conflicts:
#             for group in group:
#                 for context2 in group:
#                     if context[]:
# 
#         checklist = self.build_checklist(self.conflicts)
# 
#         for variant in self.conflicts:
#             for check in checklist[variant[2]]:
#                 if variant[0] == check[0] and variant[1] == check[1]:
#                     continue
#                 if int(variant[4]) <= int(check[5]) and int(check[4]) <= int(variant[5]):
#                     new_conflicts.append(variant)
#                     break
#             else:
#                 new_contexts.append(variant)
# =============================================================================

if __name__ == "__main__":
    import sys
    combovar = Combined_Varcon()
    combovar.add_contexts_from_file(sys.argv[1])
    combovar.divide_conflicts()
    combovar.display_counts()