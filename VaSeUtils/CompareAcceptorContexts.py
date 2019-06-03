import logging
from OverlapContext import OverlapContext

class CompareAcceptorContexts:
    def __init__(self):
        self.vaseutillogger = logging.getLogger("VaSeUtil_Logger")

    # Performs the actual comparison process
    def main(self, accconfileloc1, accconfileloc2):
        acc_contexts1 = self.read_acceptor_context_file(accconfileloc1)
        acc_contexts2 = self.read_acceptor_context_file(accconfileloc2)
        self.compare_acceptor_context_files(acc_contexts1, acc_contexts2)

    # Reads an acceptorcontext.txt file
    def read_acceptor_context_file(self, accfileloc):
        acceptor_contexts = {}
        try:
            with open(accfileloc, "r") as accfile:
                next(accfile)    # Skip the header line of the file
                for fileline in accfile:
                    fileline = fileline.strip()
                    filelinedata = fileline.split("\t")
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
            print(f"{contextid}")

        #Display info about the acceptor contexts only present in one of the two acceptor context files.
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
    def get_differences(self, difffields):
        return ""

