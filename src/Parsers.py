__author__ = 'David Magee, Carolina Correia, and Kevin Rue-Albrecht'
__copyright__ = "Copyright 2014, GPLv2"

"""Empty docstring
"""

# Module sys allows to interrup the script and return an error code
import sys
# Module gzip allows to parse compressed files without having to first decompress them.
import gzip
# Module collections allows a function to return multiple variable in a named tuple. Very clean.
import collections
# Custom module to store the data of a read
import SeqDataTypes
# Custom Module which contains a few functions for approximate matching
import ApproxMatch


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

    def next_read(self):
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

    def parse_barcode_file(self):
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
                # Do not return anything, the mapping is stored in the self.expected attribute of the BarcodesParser
                # object

    def assign_read_to_sample(self, read, barcode_length):
        """Given a read object and the expected length of a barcode, check its barcode against the expected and
        observed barcodes to define the sample it should be assigned to.

        Args:
            self, read

        Returns:
            The name of the sample to assign the read to.
        """
        # Extract the sequenced barcode from the header of the read
        sequenced_barcode = self.extract_barcode_header(read, barcode_length)
        # For training purpose, print the barcode extracted from the read header
        print("sequenced_barcode: %s\n" % sequenced_barcode)
        # if the barcode exactly matches one of the expected barcodes
        if sequenced_barcode in self.expected.keys():
            # set the assigned sample to the corresponding sample_id
            read.sample = self.expected[sequenced_barcode]
        # if the barcode exactly matches an unexpected barcodes previously resolved to a sample
        elif sequenced_barcode in self.mismatched.keys():
            # set the assigned sample to the corresponding sample_id
            read.sample = self.mismatched[sequenced_barcode]
        # if the barcode exactly matches an unexpected barcodes previously resolved to ambiguous or unmatched ones
        elif sequenced_barcode in self.unmatched:
            # set the assigned sample to False
            read.sample = False
        # if the barcode does not match any barcode observed so far
        else:
            # resolve if it can uniquely be assigned to one expected barcode
            # assigned_barcode will be the sample_id, or False if the barcode is ambiguous or unassigned
            read.sample = self.assign_to_unique_sample(sequenced_barcode)
            # If a unique sample is resolved
            if read.sample:
                # Update the observed barcodes accordingly
                self.mismatched[sequenced_barcode] = read.sample
            else:
                self.unmatched = sequenced_barcode

    def assign_to_unique_sample(self, sequenced_barcode):
        """Resolves the barcode to a unique sample, or returns False if the barcode is ambiguous.

        Args:
            self, sequenced_barcode

        Returns:
            The name of the sample matched to the barcode, or False if the barcode was ambiguous.
        """
        # Temporary list to store the expected barcodes which approximately match the sequenced barcode
        matches = []
        # For each expected barcode
        for expected in list(self.expected.keys()):
            # check if the sequenced barcode is an approximate match to each of them
            if ApproxMatch.approx_substitute(expected, sequenced_barcode, 1):
                matches.append(expected)
        # If the sequenced barcode is an approximate match of a unique expected barcode
        if len(matches) == 1:
            # return the corresponding sample name
            return matches[0]
        # If the sequenced barcode is an approximate match of more than one expected barcode, or not a match at all
        else:
            # return False
            return False

    def extract_barcode_header(self, read, barcode_length):
        """Assuming the barcode is in the header of the read between symbols # (on the left) and / (on the right)
        
        Args:
            self, read
           
        Returns:
            The sequenced barcode.
        """
        return read.header_line.split('#')[1][:barcode_length]


class AdapterParser:
    """Parses expected adapter sequence, unique identifier, and read mate to screen for it from a user file."""

    def __init__(self, filename):
        """Constructor for AdapterParser"""
        self.filename = filename

    def parse_adapters_file(self):
        """Parses the TAB-separated file and returns the adapter sequence expected in the forward mate.

        Data is expected in TAB-delimited format, with column 1 being the adapter sequence, column 2 a unique
        identifier for the adapter and column 3 one of {forward, reverse} specifying which mate to screen for this
        particular adapter sequence.

        Args:
            self

        Returns:
            An Adapter object with the adapter sequence expected in the forward read and a unique identifier.
        """
        # Initialises a list variable to temporarily store the identifier of the adapter(s) expected in the forward read
        # After parsing the file, we will check that the list contains only one element (only one type of adapter
        # should be used at each end of the RNA fragments in a given library
        forward_identifiers = []
        # Same as above for the identifiers of the adaptors expected in the reverse read
        reverse_identifiers = []
        # Initialises a variable which will contain the (last) sequence expected in the forward read
        forward_sequence = None
        # Initialises a variable which will contain the (last) sequence expected in the reverse read
        reverse_sequence = None
        # Open the file
        with open(self.filename) as file_in:
            # for each line in the file
            for line in file_in:
                # strip and split the line into a list of the different fields
                fields = line.strip().split('\t')
                # Here is the definition of the different expected fields (avoids refactoring the whole function if
                # we change our mind later)
                current_sequence = fields[0]
                current_identifier = fields[1]
                current_mate = fields[2]
                # If the current adapter is expected in the forward read
                if current_mate.lower() == "forward":
                    # add the identifier to the list of identifiers expected in the forward read
                    forward_identifiers.append(current_identifier)
                    # Store the sequence of the adapter in the prepared variable
                    forward_sequence = current_sequence
                    # We don't need to store the "mate" field: the if statement at the top of this block of code ensures
                    # that it is an adaptor expected in the forward read.
                # If the current adapter is expected in the reverse read
                elif current_mate.lower() == "reverse":
                    # add the identifier to the list of identifiers expected in the reverse read
                    reverse_identifiers.append(current_identifier)
                    # Store the sequence of the adapter in the prepared variable
                    reverse_sequence = current_sequence
                    # We don't need to store the "mate" field: the if statement at the top of this block of code ensures
                    # that it is an adaptor expected in the forward read.
                # If the current adapter has an unexpected value in the "mate" field
                else:
                    # return an error and exit the script
                    print(
                        'Unexpected value in the "mate" field. Please check your adaptor file (%s) or seek advice '
                        'about this script!' % self.filename)
                    sys.exit(6)
        # Once we finished parsing the adaptor file, we do some sanity checks and if ok, return the adaptor object
        # return an error and exit the script if there are more than 1 adaptor expected in the forward read
        if len(forward_identifiers) > 1:
            print(
                'More than one adaptor sequence is expected in the forward read. Please check your adaptor file (%s) '
                'or seek advice about this script!' % self.filename)
            sys.exit(7)
            # Return an error
        if len(reverse_identifiers) > 1:
            print(
                'More than one adaptor sequence is expected in the reverse read. Please check your adaptor file (%s) '
                'or seek advice about this script!' % self.filename)
            sys.exit(8)
        # Return an error if the identifier from forward and reverse read are the same (if we reached this line,
        # we know there is one unique adapter for each mate).
        if forward_identifiers == reverse_identifiers:
            print(
                'Identifiers from forward and reverse are not unique. Please check your adaptor file (%s) or seek '
                'advice about this script!' % self.filename)
            sys.exit(9)
        # After all the sanity check return the forward and reverse adaptor objects as a tuple
        # Prepare a named tuple structure to cleanly return the two adapters in one variable
        adapter_pair = collections.namedtuple('AdapterPair', ['forward', 'reverse'])
        # Fill in the tuple structure with the adapter objects and return it
        return adapter_pair(forward=SeqDataTypes.Adapter(forward_identifiers[0], forward_sequence, 'forward'),
                            reverse=SeqDataTypes.Adapter(reverse_identifiers[0], reverse_sequence, 'reverse'))