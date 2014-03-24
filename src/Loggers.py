__author__ = 'David Magee, Carolina Correia, and Kevin Rue-Albrecht'
__copyright__ = "Copyright 2014, GPLv2"

"""Empty docstring
"""

# Module re allows to match regular expressions. Useful to see if a filename comes from BIG, MSU, or Conway.
import re


class ReadLogger:
    """Tracks number of reads assigned to each sample, reads filtered out because of quality or adapter
    contamination, etc."""

    def __init__(self, filename, samples):
        """Constructor for ReadLogger
        
        Note that filename is the name of the input forward read file, not the name of the report file. The latter 
        will be dynamically created when initialising the logger.
        """
        # Name pf the file to write the statistics to
        self.filename = set_report_filename_from_raw_file(filename)
        print("Info: Report file: %s" % self.filename)
        # Total number of read pairs processed
        # self.read_count = 0
        # = the sum of all assigned + unassigned reads. No need to do an extra counter for that!
        # Number of reads assigned to each barcode/sample
        self.assigned = dict.fromkeys(samples, 0)
        #print("Test: assigned: ", self.assigned)
        # Number of unassigned reads
        self.unassigned = 0
        # Number of reads excluded because of adapter presence per sample
        self.adapter_excluded = dict.fromkeys(samples, 0)
        # Number of reads discarded because of poor quality per sample
        self.quality_excluded = dict.fromkeys(samples, 0)
        # Number of reads passed filtering per sample
        self.accepted = dict.fromkeys(samples, 0)


def set_report_filename_from_raw_file(filename):
    """Generates a report filename based on the name of the filename of the raw reads.

    Note that the filename of the raw reads will look different depending on the sequencing centre. For now only the
    BGI case will be implemented to test the AlvMac files.

    Args:
        filename

    Returns:
        None. Sets the self.filename field instead.
    """
    # For each know centre, test if the filename matches its typical pattern
    # Test for BGI (Raw_"anything"_1."anything")
    # keep the "Raw_anything" until the character before the _1
    BGI_match = re.match('Raw_.*(?=_1\..*$)', filename)
    test_match = re.match('single(?=_pe1*)|tenthousands(?=_pe1*)', filename)
    if BGI_match:
        print("Info: Sequencing centre: BGI (Based on raw reads filename)")
        report_basename = BGI_match.group(0)
    # Test for our training files
    # keep the name from the start until the character before the _pe1
    elif test_match:
        print("Info: Sequencing centre: Testing data (Based on raw reads filename)")
        report_basename = test_match.group(0)
    # if no pattern matches
    else:
        # return an error and exit
        print(
            "Error: Your data seems to have been generated by an unknown sequencing centre; please seek advice about "
            "this script!")
        sys.exit(2)
    # if we reach here, a basename was defined, complete it with the extension
    return report_basename + "_1_2.report"