import logging
from OverlapContext import OverlapContext

class CompareAcceptorContexts:
    def __init__(self, vaseutilhelper):
        self.vaseutillogger = logging.getLogger("VaSeUtil_Logger")
        self.vuh = vaseutilhelper

    # Performs the actual comparison process
    def main(self, accconfileloc1, accconfileloc2):
        acc_contexts1 = self.read_acceptor_context_file(accconfileloc1)
        acc_contexts2 = self.read_acceptor_context_file(accconfileloc2)
        self.compare_acceptor_context_files(acc_contexts1, acc_contexts2)

    # Reads an acceptorcontext.txt file
    def read_acceptor_context_file(self, accfileloc, samplefilter=None, contextfilter=None, chromfilter=None):
        acceptor_contexts = {}
        try:
            with open(accfileloc, "r") as accfile:
                next(accfile)    # Skip the header line of the file
                for fileline in accfile:
                    fileline = fileline.strip()
                    filelinedata = fileline.split("\t")

                    # Check whether the acceptor context passes any set filters
                    samplepass = self.vuh.passes_filter(filelinedata[1], samplefilter)
                    contextpass = self.vuh.passes_filter(filelinedata[0], contextfilter)
                    chrompass = self.vuh.passes_filter(filelinedata[2], chromfilter)

                    if samplepass and contextpass and chrompass:
                        acc_context = OverlapContext(filelinedata[0], filelinedata[1], filelinedata[2], filelinedata[3],
                                                     filelinedata[4], filelinedata[5], filelinedata[7].split(";"))
                        acceptor_contexts[filelinedata[0]] = acc_context
        except IOError:
            self.vaseutillogger.warning(f"Could not read acceptor context file: {accfileloc}")
        return acceptor_contexts

    # Performs the comparison of the two acceptor context files
    def compare_acceptor_context_files(self, acc_contexts1, acc_contexts2):
        shared_contexts = list(set(acc_contexts1.keys()) & set(acc_contexts2.keys()))
        for contextid in shared_contexts:
            varcon_compare_results = acc_contexts1[contextid].compare(acc_contexts2[contextid])
            if len(varcon_compare_results) > 0:
                print(f"{contextid}: {self.get_differences(varcon_compare_results)}")

        # Display info about the acceptor contexts only present in one of the two acceptor context files.
        self.unique_acceptor_context_info(acc_contexts1, acc_contexts2)
        self.unique_acceptor_context_info(acc_contexts2, acc_contexts1)

    # Displays information about acceptor contexts only present in one of the two acceptor context files
    def unique_acceptor_context_info(self, acc_contexts1, acc_contexts2):
        unique_contexts = list(set(acc_contexts1.keys()) - set(acc_contexts2.keys()))
        for contextid in unique_contexts:
            diffs = acc_contexts1[contextid].compare(acc_contexts2[contextid])
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

