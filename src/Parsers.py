__author__ = 'David Magee, Carolina Correia, and Kevin Rue-Albrecht'
__copyright__ = "Copyright 2014, GPLv2"

"""Empty docstring
"""

# Module sys allows to interrup the script and return an error code
import sys
# Module gzip allows to parse compressed files without having to first decompress them.
import gzip
# Custom module to store the data of a read
import SeqDataTypes


class FastqgzParser:
    """Parses compressed fastq files."""

    def __init__(self, filename):
        """Constructor for FastqqzParser
        :rtype : FastqgzParser
        """
        self.filename = filename
        self.file_handler = gzip.open(self.filename, 'rt')

    def open(self):
        """Opens the file stream.

        Args:
            self

        Returns:
            None
        """
        self.file_handler = gzip.open(self.filename, "rt")

    def nextRead(self):
        """Parses and return the next sequenced read in the file.

        Args:
            self

        Returns:
            The next sequenced read in the file.
        """
        header_line = self.file_handler.readline().strip()
        sequence_line = self.file_handler.readline().strip()
        separator_line = self.file_handler.readline().strip()
        quality_line = self.file_handler.readline().strip()
        return SeqDataTypes.Read(header_line, sequence_line, separator_line, quality_line)

    # The first listed Read is a module. The second listed read is an objectClass type Read. This code takes each of
    # the four variables (header_line, sequence_line, separator_line, quality_line) and inserts them into one single
    # variable.

    def close(self):
        """Closes the file stream.
        
        Args:
            self
           
        Returns:
            None
        """
        self.file_handler.close()


class BarcodesParser:
    """Parses expected barcodes in a file and observed barcodes in the reads."""

    def __init__(self, filename):
        """Constructor for BarcodesParser"""
        self.filename = filename
        self.expected = {}
        self.mismatched = {}
        self.unmatched = []

    def parseBarcodeFile(self):
        """Parses the file provided with the expected barcodes and saves them in the self.expected variable.

        The samples_identifiers (sample_ids) and barcodes will first be stored in temporary lists to check for issues
        like duplicates (all identifiers and barcodes need to be unique). The script will be interrupted if issues are
        found. Otherwise, the pairs will be saved in a dictionary.

        Args:
            self

        Returns:
            None.
        """
        # Initialise a list variable to temporarily store the sample identifiers
        sample_ids = []
        # Initialise a list variable to temporarily store the barcodes
        barcodes = []
        # Open the user file for import
        with open(self.filename) as file_in:
            # For each line in the file
            for line in file_in:
                # strip and split the line into the tab-separated fields
                fields = line.strip().split('\t')
                # Append the current barcode to the list of barcodes
                barcodes.append(fields[0])
                # Append the current sample_id to the list of sample_ids
                sample_ids.append(fields[1])
                # Note that the sample_ids and barcodes are stored in the same order, therefore the first brcode
                # corresponds to the first sample_id, and so on
            # After all the lines of the file are read, do the sanity checks, and then store the mapping in a dictionary
            # Return an error if there are duplicate sample_ids (who puts two barcodes for the same sample?)
            if len(set([x for x in sample_ids if sample_ids.count(x) > 1])):
                print(
                    'The same sample identifier appears twice in the barcode file. Please check again or seek advice'
                    ' about this script!\n')
                sys.exit(5)
            # Return an error if there are duplicate barcodes (who puts two samples with the same barcode?)
            if len(set([x for x in barcodes if barcodes.count(x) > 1])):
                print(
                    'The same sample barcode appears twice in the barcode file. Please check again or seek advice'
                    ' about this script!\n')
                sys.exit(6)
            # When all the checks are done, store all the barcode->sample pairs in a dictionary for deconvolution
            for index in range(len(barcodes)):
                self.expected[barcodes[index]] = sample_ids[index]
        # Do not return anything, the mapping is stored in the self.expected attribute of the BarcodesParser object