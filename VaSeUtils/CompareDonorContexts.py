import logging
from OverlapContext import OverlapContext

class CompareDonorContexts:
    def __init__(self, vaseutilhelper):
        self.vaseutillogger = logging.getLogger("VaSeUtil_Logger")
        self.vuh = vaseutilhelper

    def main(self, donconfileloc1, donconfileloc2):
        don_contexts1 = self.read_donor_context_file(donconfileloc1)
        don_contexts2 = self.read_donor_context_file(donconfileloc2)

        self.compare_donor_context_files(don_contexts1, don_contexts2)

    # Reads the donor context file
    def read_donor_context_file(self, contextfileloc):
        donor_contexts = {}
        try:
            with open(contextfileloc, 'r') as dcontextfile:
                next(dcontextfile)    # Skip the header line.
                for fileline in dcontextfile:
                    fileline = fileline.strip()
                    filelinedata = fileline.split("\t")
                    don_context = OverlapContext(filelinedata[0], filelinedata[1], filelinedata[2], filelinedata[3],
                                                 filelinedata[4], filelinedata[5], filelinedata[7].split(";"))
                    donor_contexts[filelinedata[0]] = don_context

        except IOError as ioe:
            self.vaseutillogger.warning("Could not read donor context file " +str(ioe.filename))
        return donor_contexts

    # Compares two donor context files
    def compare_donor_context_files(self, don_contexts1, don_contexts2):
        shared_contexts = list(set(don_contexts1.keys()) & set(don_contexts2.keys()))
        for contextid in shared_contexts:
            varcon_compare_results = don_contexts1[contextid].compare(don_contexts2[contextid])
            if len(varcon_compare_results) > 0:
                print(f"{contextid}: {self.get_differences(varcon_compare_results)}")

        # Display info about the acceptor contexts only present in one of the two acceptor context files.
        self.unique_donor_context_info(don_contexts1, don_contexts2)
        self.unique_donor_context_info(don_contexts2, don_contexts1)

    # Displays information about acceptor contexts only present in one of the two acceptor context files
    def unique_donor_context_info(self, don_contexts1, don_contexts2):
        unique_contexts = list(set(don_contexts1.keys()) - set(don_contexts2.keys()))
        for contextid in unique_contexts:
            diffs = don_contexts1[contextid].compare(don_contexts2[contextid])
            if len(diffs) > 0:
                print(self.get_differences(diffs))

    # Displays the differences between two shared acceptor contexts
    def get_differences(self, diffields):
        difmap = self.vuh.get_accdon_context_fields()
        differences = ""
        for diffield in diffields:
            if diffield in difmap:
                differences += f"{difmap[diffield]}: {diffields[diffield][0]}/{diffields[diffield][1]}"
        return differences
