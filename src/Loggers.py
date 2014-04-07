__author__ = 'David Magee, Carolina Correia, and Kevin Rue-Albrecht'
__copyright__ = "Copyright 2014, GPLv2"

"""Empty docstring
"""

# Module os.path allows to extract the basename of a file
import os.path
# Module re allows to match regular expressions. Useful to see if a filename comes from BIG, MSU, or Conway.
import re
# Module sys allows to interrupt the script with an error code.
import sys

class ReadLogger:
    """Tracks number of reads assigned to each sample, reads filtered out because of quality or adapter
    contamination, etc.

    This type of object was given variables and functions to count the number of reads with an ambiguous barcode,
    adapter contamination, poor quality, or successfully passing all quality control steps."""

    def __init__(self, filename, barcodes_parser, basename_raw):
        """Constructor for ReadLogger
        
        Note that filename is the name of the input forward read file, not the name of the report file. The latter 
        will be dynamically created when initialising the logger.

        Note that barcodes_samples is a dictionary mapping each expected barcode (key) to the corresponding sample (
        value).
        """
        # Name of the file to write_stats the statistics to
        self.filename = basename_raw + "_1_2.report.txt"
        print("Info: Report file: %s" % self.filename)
        # Saves the link to the BarcodesParser, it contains the mapping between barcodes and samples, useful to write
        #  the report file
        self.barcode_parser = barcodes_parser
        # Total number of read pairs processed
        # self.read_count = 0
        # = the sum of all assigned + unassigned reads. No need to do an extra counter for that!
        # Number of reads assigned to each barcode/sample
        self.assigned = dict.fromkeys(barcodes_parser.expected.values(), 0)
        #print("Test: assigned: ", self.assigned)
        # Number of unassigned reads
        self.unassigned = 0
        # Number of reads excluded because of adapter presence per sample
        self.adapter_excluded = dict.fromkeys(barcodes_parser.expected.values(), 0)
        # Number of reads discarded because of poor quality per sample
        self.quality_excluded = dict.fromkeys(barcodes_parser.expected.values(), 0)
        # Number of reads passed filtering per sample
        self.accepted = dict.fromkeys(barcodes_parser.expected.values(), 0)

    def write_stats(self):
        """Writes the current statistics in the report file.

        Args:
            self

        Returns:
            None
        """
        # Opens the report file
        #print("Test: expected barcodes:", self.barcode_parser.expected)
        with open(self.filename, "w") as report:
            # write the header line
            report.write(
                "Barcode\tSample\tPre-filtering\tPost-filtering\t% post-filtering\tAdapter filtered\t% adapter "
                "filtered\tQuality filtered\t% quality filtered\n")
            # for each expected barcode
            for barcode in self.barcode_parser.expected.keys():
                #print("Test: expected_barcodes:", barcode.decode("ascii"))
                # get the corresponding sample name (instead of getting it everytime below)
                sample_name = self.barcode_parser.expected[barcode]
                # Pre-calculates the ratios for clarity and to deal with the potential division by 0
                ratio_accepted = self.accepted[sample_name] / self.assigned[sample_name] if self.assigned[
                    sample_name]  else 0
                ratio_adapter = self.adapter_excluded[sample_name] / self.assigned[sample_name] if self.assigned[
                    sample_name]  else 0
                ratio_quality = self.quality_excluded[sample_name] / self.assigned[sample_name] if self.assigned[
                    sample_name]  else 0
                # write the statistics
                report.write("%s\t%s\t%i\t%i\t%.3f\t%i\t%.3f\t%i\t%.3f\n" % (barcode.decode("ascii"),
                                                                             sample_name,
                                                                             self.assigned[sample_name],
                                                                             self.accepted[sample_name],
                                                                             ratio_accepted,
                                                                             self.adapter_excluded[sample_name],
                                                                             ratio_adapter,
                                                                             self.quality_excluded[sample_name],
                                                                             ratio_quality))
            # add a separator line to separate from merged statistics
            report.write("\n")
            ratio_accepted = sum(self.accepted.values()) / sum(self.assigned.values()) if sum(
                self.assigned.values()) else 0
            ratio_adapter = sum(self.adapter_excluded.values()) / sum(self.assigned.values()) if sum(
                self.assigned.values()) else 0
            ratio_quality = sum(self.quality_excluded.values()) / sum(self.assigned.values()) if sum(
                self.assigned.values()) else 0
            # add a line for the assigned total statistics
            report.write("Assigned\t\t%i\t%i\t%.3f\t%i\t%.3f\t%i\t%.3f\n" % (sum(self.assigned.values()),
                                                                             sum(self.accepted.values()),
                                                                             ratio_accepted,
                                                                             sum(self.adapter_excluded.values()),
                                                                             ratio_adapter,
                                                                             sum(self.quality_excluded.values()),
                                                                             ratio_quality))
            # add a line for the unassigned reads statistics
            report.write("Unassigned\t\t%i\n" % self.unassigned)
            # add aline for the total read
            report.write("Total\t\t%i\n" % (sum(self.assigned.values()) + self.unassigned))
