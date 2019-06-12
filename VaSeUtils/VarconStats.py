import logging


class VarconStats:
    def __init__(self, varconid, avg_arlength, avg_drlength, med_arlength, med_drlength, avg_arqual, avg_drqual,
                 med_arqual, med_drqual, avg_amapq, avg_dmapq, med_amapq, med_dmapq):
        self.context_id = varconid
        self.avg_aread_length = avg_arlength
        self.avg_dread_length = avg_drlength
        self.median_aread_length = med_arlength
        self.median_dread_length = med_drlength
        self.avg_aread_qual = avg_arqual
        self.avg_dread_qual = avg_drqual
        self.median_aread_qual = med_arqual
        self.median_dread_qual = med_drqual
        self.avg_amapq = avg_amapq
        self.avg_dmapq = avg_dmapq
        self.median_amapq = med_amapq
        self.median_dmapq = med_dmapq
        self.vaseutillogger = logging.getLogger("VaSeUtil_Logger")

    # Returns the variant context identifier
    def get_variant_context_id(self):
        return self.context_id

    # Returns the average acceptor read length of the variant context
    def get_average_acceptor_read_length(self):
        return self.avg_aread_length

    # Returns the average donor read length of the variant context
    def get_average_donor_read_length(self):
        return self.avg_dread_length

    # Returns the acceptor/donor average read length ratio
    def get_average_read_length_ratio(self):
        return float(self.avg_aread_length / self.avg_dread_length)

    # Returns the difference between the average acceptor and donor length
    def get_average_read_length_difference(self):
        return float(self.avg_aread_length - self.avg_dread_length)

    # Returns the median acceptor read length of the variant context
    def get_median_acceptor_read_length(self):
        return self.median_aread_length

    # Returns the median donor read length of the variant context
    def get_median_donor_read_length(self):
        return self.median_dread_length

    # Returns the acceptor/donor median read length ratio
    def get_median_read_length_ratio(self):
        return float(self.median_aread_length / self.median_dread_length)

    # Returns the difference between the acceptor/donor median read length
    def get_median_read_length_difference(self):
        return float(self.median_aread_length - self.median_dread_length)

    # Returns the average acceptor read quality of the variant context
    def get_average_acceptor_read_quality(self):
        return self.avg_aread_qual

    # Returns the average donor read quality of the variant context
    def get_average_donor_read_quality(self):
        return self.avg_dread_qual

    # Returns the acceptor/donor average read quality
    def get_average_read_quality_ratio(self):
        return float(self.avg_aread_qual / self.avg_dread_qual)

    # Returns the acceptor/donor average read quality
    def get_average_read_quality_difference(self):
        return float(self.avg_aread_qual - self.avg_dread_qual)

    # Returns the median acceptor read quality of the variant context
    def get_median_acceptor_read_quality(self):
        return self.median_aread_qual

    # Returns the median donor read quality of the variant context
    def get_median_donor_read_quality(self):
        return self.median_dread_qual

    # Returns the average acceptor mapq values of the variant context
    def get_average_acceptor_mapq(self):
        return self.avg_amapq

    # Returns the average donor mapq values of the variant context
    def get_average_donor_mapq(self):
        return self.avg_dmapq

    # Returns the median acceptor mapq of the variant context
    def get_median_acceptor_mapq(self):
        return self.median_amapq

    # Returns the median donor mapq of the variant context
    def get_median_donor_mapq(self):
        return self.median_dmapq

    # Returns the variant context statistics as a file entry
    def to_string(self):
        return str(self.context_id) + "\t" + str(varconAvgAReadLength) + "\t" + str(varconAvgDReadLength) + "\t" + str(varconMedianAReadLength) + "\t" + str(varconMedianDReadLength) + "" + str() + "" + str() + "" + str() + "" + str() + "" + str() + ""
