__author__ = 'David Magee, Carolina Correia, and Kevin Rue-Albrecht'
__copyright__ = "Copyright 2014, GPLv2"

"""Empty docstring
"""

# Module re allows to match regular expressions. Useful to see if a filename comes from BIG, MSU, or Conway.
import re


class ReadLogger:
    """Tracks number of reads assigned to each sample, reads filtered out because of quality or adapter
    contamination, etc."""

    def __init__(self, filename, barcodes_parser):
        """Constructor for ReadLogger
        
        Note that filename is the name of the input forward read file, not the name of the report file. The latter 
        will be dynamically created when initialising the logger.

        Note that barcodes_samples is a dictionary mapping each expected barcode (key) to the corresponding sample (
        value).
        """
        # Name of the file to write_stats the statistics to
        self.filename = set_report_filename_from_raw_file(filename)
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
        with open(self.filename, "w") as report:
            # write the header line
            report.write(
                "Barcode\tSample\tPre-filtering\tPost-filtering\t% post-filtering\tAdapter filtered\t% adapter "
                "filtered\tQuality filtered\t% quality filtered\n")
            # for each expected barcode
            for barcode in self.barcode_parser.expected.keys():
                # get the corresponding sample name (instead of getting it everytime below)
                sample_name = self.barcode_parser.expected[barcode]
                # Pre-calculates the ratios for clarity and to deal with the potential division by 0
                ratio_accepted = self.accepted[sample_name] / self.assigned[sample_name] if self.assigned[
                    sample_name]  else 0
                ratio_adapter = self.quality_excluded[sample_name] / self.assigned[sample_name] if self.assigned[
                    sample_name]  else 0
                ratio_quality = self.adapter_excluded[sample_name] / self.assigned[sample_name] if self.assigned[
                    sample_name]  else 0
                # write the statistics
                report.write("%s\t%s\t%i\t%i\t%.3f\t%i\t%.3f\t%i\t%.3f\n" % (barcode,
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
            report.write("Total\t%i\n" % (sum(self.assigned.values()) + self.unassigned))


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
    return report_basename + "_1_2.report.txt"
