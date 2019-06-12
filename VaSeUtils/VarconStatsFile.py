import logging
from VarconStats import VarconStats


class VarconStatsFile:
    def __init__(self, statsfileloc, contextfilter=None):
        self.vaseutillogger = logging.getLogger("VaSeUtil_Logger")
        self.varcon_stats_data = self.read_varcon_stats_file(statsfileloc, contextfilter)

    # Reads the variant context statistics file.
    def read_varcon_stats_file(self, fileloc, contextfilter=None):
        varcon_stats = {}
        try:
            with open(fileloc, 'r') as varconstatsfile:
                next(varconstatsfile)    # Skips the header line
                for fileline in varconstatsfile:
                    fileline = fileline.strip()
                    filelinedata = fileline.split("\t")
                    if self.passes_filter(filelinedata[0], contextfilter):
                        varcon_stats[filelinedata[0]] = VarconStats(filelinedata[0], filelinedata[1], filelinedata[2],
                                                                    filelinedata[3], filelinedata[4], filelinedata[5],
                                                                    filelinedata[6], filelinedata[7], filelinedata[8],
                                                                    filelinedata[9], filelinedata[10], filelinedata[11],
                                                                    filelinedata[12])
        except IOError as ioe:
            self.vaseutillogger.warning(f"Could not read varconstats file {ioe.filename}")
        return varcon_stats

    # Writes the variant context statistics file to a specified output file
    def write_varcon_stats_file(self, outfileloc, contextfilter=None):
        try:
            with open(outfileloc, 'w') as outfile:
                for contextid in self.varcon_stats_data:
                    if self.passes_filter(contextid, contextfilter):
                        outfile.write(f"{self.varcon_stats_data[contextid].to_string()}\n")
        except IOError:
            self.vaseutillogger.warning(f"Could not write varconstats file to {outfileloc}")

    # Returns the entire variant context statistics data map
    def get_varcon_stats_data(self):
        return self.varcon_stats_data

    # Returns an entire VarconStats object containing all statistics for a single variant context id.
    def get_varcon_stats(self, contextid):
        if contextid in self.varcon_stats_data:
            return 
        return None

    # Returns the variant context identifier
    def get_variant_context_id(self, contextid):
        if contextid in self.varcon_stats_data:
            return self.varcon_stats_data[contextid].get_variant_context_id()
        return None

    # Returns the average acceptor read length of the variant context
    def get_average_acceptor_read_length(self, contextid):
        if contextid in self.varcon_stats_data:
            return self.varcon_stats_data[contextid].get_average_acceptor_read_length()
        return None

    # Returns the average donor read length of the variant context
    def get_average_donor_read_length(self, contextid):
        if contextid in self.varcon_stats_data:
            return self.varcon_stats_data[contextid].get_average_donor_read_length()
        return None

    # Returns the median acceptor read length of the variant context
    def get_median_acceptor_read_length(self, contextid):
        if contextid in self.varcon_stats_data:
            return self.varcon_stats_data[contextid].get_median_acceptor_read_length()
        return None

    # Returns the median donor read length of the variant context
    def get_median_donor_read_length(self, contextid):
        if contextid in self.varcon_stats_data:
            return self.varcon_stats_data[contextid].get_median_donor_read_length()
        return None

    # Returns the average acceptor read quality of the variant context
    def get_average_acceptor_read_quality(self, contextid):
        if contextid in self.varcon_stats_data:
            return self.varcon_stats_data[contextid].get_average_acceptor_read_quality()
        return None

    # Returns the average donor read quality of the variant context
    def get_average_donor_read_quality(self, contextid):
        if contextid in self.varcon_stats_data:
            return self.varcon_stats_data[contextid].get_average_donor_read_quality()
        return None

    # Returns the median acceptor read quality of the variant context
    def get_median_acceptor_read_quality(self, contextid):
        if contextid in self.varcon_stats_data:
            return self.varcon_stats_data[contextid].get_median_acceptor_read_quality()
        return None

    # Returns the median donor read quality of the variant context
    def get_median_donor_read_quality(self, contextid):
        if contextid in self.varcon_stats_data:
            return self.varcon_stats_data[contextid].get_median_donor_read_quality()
        return None

    # Returns the average acceptor mapq values of the variant context
    def get_average_acceptor_mapq(self, contextid):
        if contextid in self.varcon_stats_data:
            return self.varcon_stats_data[contextid].get_average_acceptor_mapq()
        return None

    # Returns the average donor mapq values of the variant context
    def get_average_donor_mapq(self, contextid):
        if contextid in self.varcon_stats_data:
            return self.varcon_stats_data[contextid].get_average_donor_mapq()
        return None

    # Returns the median acceptor mapq of the variant context
    def get_median_acceptor_mapq(self, contextid):
        if contextid in self.varcon_stats_data:
            return self.varcon_stats_data[contextid].get_median_acceptor_mapq()
        return None

    # Returns the median donor mapq of the variant context
    def get_median_donor_mapq(self, contextid):
        if contextid in self.varcon_stats_data:
            return self.varcon_stats_data[contextid].get_median_donor_mapq()
        return None

    # Returns a list of all variant contexts
    def get_variant_contexts_ids(self):
        return list(self.varcon_stats_data.keys())
    
    # Returns a list/hashmap of all average acceptor read lengths
    def get_average_acceptor_read_lengths(self, aslist=False):
        if aslist:
            return [x.get_average_acceptor_read_length() for x in list(self.varcon_stats_data.values())]
        return {k: v.get_average_acceptor_read_length() for k, v in self.varcon_stats_data.items()}

    # Returns a list/hashmap of all average donor read lengths
    def get_average_donor_read_lengths(self, aslist=False):
        if aslist:
            return [x.get_average_donor_read_length() for x in list(self.varcon_stats_data.values())]
        return {k: v.get_average_donor_read_length() for k, v in self.varcon_stats_data.items()}

    # Returns a list/hashmap of all median ascceptor read lengths
    def get_median_acceptor_read_lengths(self, aslist=False):
        if aslist:
            return [x.get_median_acceptor_read_length() for x in list(self.varcon_stats_data.values())]
        return {k: v.get_median_acceptor_read_length() for k, v in self.varcon_stats_data.items()}

    # Returns a list/hashmap of all median donor read lengths
    def get_median_donor_read_lengths(self, aslist=False):
        if aslist:
            return [x.get_median_donor_read_length() for x in list(self.varcon_stats_data.values())]
        return {k: v.get_median_donor_read_length() for k, v in self.varcon_stats_data.items()}

    # Returns a list/hashmap of all average acceptor read qualities
    def get_average_acceptor_read_qualities(self, aslist=False):
        if aslist:
            return [x.get_average_acceptor_read_quality() for x in list(self.varcon_stats_data.values())]
        return {k: v.get_average_acceptor_read_quality() for k, v in self.varcon_stats_data.items()}

    # Returns a list/hashmap of all average donor read qualities
    def get_average_donor_read_qualities(self, aslist=False):
        if aslist:
            return [x.get_average_donor_read_quality() for x in list(self.varcon_stats_data.values())]
        return {k: v.get_average_donor_read_quality() for k, v in self.varcon_stats_data.items()}

    # Returns a list/hashmap of all median acceptor read qualities
    def get_median_acceptor_read_qualities(self, aslist=False):
        if aslist:
            return [x.get_median_acceptor_read_quality() for x in list(self.varcon_stats_data.values())]
        return {k: v.get_median_acceptor_read_quality() for k, v in self.varcon_stats_data.items()}

    # Returns a list/hashmap of all median donor read qualities
    def get_median_donor_read_qualities(self, aslist=False):
        if aslist:
            return [x.get_median_donor_read_quality() for x in list(self.varcon_stats_data.values())]
        return {k: v.get_median_donor_read_quality() for k, v in self.varcon_stats_data.items()}

    # Returns a list/hashmap of all average acceptor read mapq values
    def get_average_acceptor_mapqs(self, aslist=False):
        if aslist:
            return [x.get_average_acceptor_mapq() for x in list(self.varcon_stats_data.values())]
        return {k: v.get_average_acceptor_mapq() for k, v in self.varcon_stats_data.items()}

    # Returns a list/hashmap of all average donor read mapq values
    def get_average_donor_mapqs(self, aslist=False):
        if aslist:
            return [x.get_average_donor_mapq() for x in list(self.varcon_stats_data.values())]
        return {k: v.get_average_donor_mapq() for k, v in self.varcon_stats_data.items()}

    # Returns a list/hashmap of all median acceptor read mapq values
    def get_median_acceptor_mapqs(self, aslist=False):
        if aslist:
            return [x.get_median_acceptor_mapq() for x in list(self.varcon_stats_data.values())]
        return {k: v.get_median_acceptor_mapq() for k, v in self.varcon_stats_data.items()}
    
    # Returns a list/hashmap of all median donor read mapq values
    def get_median_donor_mapqs(self, aslist=False):
        if aslist:
            return [x.get_median_donor_mapq() for x in list(self.varcon_stats_data.values())]
        return {k: v.get_median_donor_mapq() for k, v in self.varcon_stats_data.items()}

    # Checks whether a value is within the filter
    def passes_filter(self, valtocheck, filterlist):
        if filterlist is not None:
            return valtocheck in filterlist
        return True
