__author__ = 'David Magee, Carolina Correia, and Kevin Rue-Albrecht'
__copyright__ = "Copyright 2014, GPLv2"

"""Empty docstring
"""

# Module collections allows a function to return multiple variable in a named tuple. Very clean.
import collections
# Module fuzzysearch allows approximate matching of a substring within another string
import fuzzysearch
# Module gzip allows to parse compressed files without having to first decompress them.
import gzip
# Module re allows to match regular expressions. Useful to see if a filename comes from BIG, MSU, or Conway.
import re
# Module os.path allows to extract the basename of a file
import os.path
# Module sys allows to interrupt the script and return an error code
import sys
# Custom module to store the data of a read
import SeqDataTypes
# Custom Module which contains a few functions for approximate matching
import ApproxMatch


class FastqgzPairParser:
    """Synchronises the parsing of a forward and a reverse compressed fastq files."""

    def __init__(self, forward_file, reverse_file, expected_barcode_pattern):
        """Constructor for ReadPairParser"""
        self.forward_parser = FastqgzParser(forward_file, expected_barcode_pattern)
        self.reverse_parser = FastqgzParser(reverse_file, expected_barcode_pattern)

    def next_read_pair(self):
        """Reads the next mate pair.

        Args:
            self

        Returns:
            An object of class ReadPair containing the two reads as objects themselves.
        """
        # Fill in the tuple structure with the adapter objects and return it
        return SeqDataTypes.ReadPair(self.forward_parser.next_read(), self.reverse_parser.next_read())


class FastqgzParser:
    """Parses compressed fastq files."""

    def __init__(self, filename, expected_barcode_pattern):
        """Constructor for FastqqzParser
        :rtype : FastqgzParser
        """
        self.filename = filename
        self.barcode_pattern = expected_barcode_pattern
        self.open()

    def open(self):
        """Opens the file stream.

        Args:
            self

        Returns:
            None
        """
        lines = [
            line.strip() for line in
            gzip.open(self.filename, "rb").read().splitlines()
        ]
        def next_read_iter(lines=lines):
            for index in range(0, len(lines), 4):
                yield SeqDataTypes.Read(lines[index], lines[index + 1], lines[index + 2], lines[index + 3])
            while True:
                yield SeqDataTypes.Read(b'', b'', b'', b'')
        self._next_read_iter = next_read_iter()

    def next_read(self):
        """Parses and return the next sequenced read in the file.

        Args:
            self

        Returns:
            The next sequenced read in the file.
        """
        # Will crash if the expected barcodes are longer than the sequence between # and /
        # Return a Read object
        return next(self._next_read_iter)

    # The first listed Read is a module. The second listed read is an objectClass type Read. This code takes each of
    # the four variables (header_line, sequence_line, separator_line, quality_line) and inserts them into one single
    # variable.


class BarcodesParser:
    """Parses expected barcodes in a file and observed barcodes in the reads.

    # This type of object was given variables and functions to import expected barcodes from a user file. It also
    creates and updates a list of barcodes observed in the processed reads, each of which is resolved once and
    for all to either a unique sample or declared ambiguous and leads to the immediate exclusion of all reads
    showing that same barcode
    ."""

    def __init__(self, filename):
        """Constructor for BarcodesParser"""
        self.filename = filename
        self.expected = {}
        self.mismatched = {}
        self.unmatched = []
        self.barcode_length = None
        self.barcode_pattern = None

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
                barcodes.append(fields[0].encode('ascii'))
                # Append the current sample_id to the list of sample_ids
                sample_ids.append(fields[1])
                # Note that the sample_ids and barcodes are stored in the same order, therefore the first brcode
                # corresponds to the first sample_id, and so on
        # After all the lines of the file are read, do the sanity checks, and then store the mapping in a dictionary
        # Return an error if there are duplicate sample_ids (who puts two barcodes for the same sample?)
        if len(set([x for x in sample_ids if sample_ids.count(x) > 1])):
            print(
                'Error: The same sample identifier appears twice in the barcode file. Please check again or seek '
                'advice about this script!\n')
            sys.exit(5)
        # Return an error if there are duplicate barcodes (who puts two samples with the same barcode?)
        if len(set([x for x in barcodes if barcodes.count(x) > 1])):
            print(
                'Error: The same sample barcode appears twice in the barcode file. Please check again or seek '
                'advice about this script!\n')
            sys.exit(6)
        # When all the checks are done, store all the barcode->sample pairs in a dictionary for deconvolution
        for index in range(len(barcodes)):
            self.expected[barcodes[index]] = sample_ids[index]
            # Do not return anything, the mapping is stored in the self.expected attribute of the BarcodesParser
            # object
        # For training purpose, print the dictionary of expected barcodes
        #print("Test: expected_barcodes:", self.expected)
        # From the first barcode in the variable, define the expected barcode length, assuming all barcodes are the same
        # length
        self.barcode_length = len(list(self.expected.keys())[0])
        # Assuming the barcode is in the header after # symbol, compile once and for all a motif/pattern
        self.barcode_pattern = re.compile(("#([ATGCN]{%i})[ATGCN]*/" % self.barcode_length).encode('ascii'))

    
    def assign_read_to_sample(self, read_pair, barcode_length):
        """Given a read object and the expected length of a barcode, check its barcode against the expected and
        observed barcodes to define the sample it should be assigned to.

        Args:
            self, read

        Returns:
            The name of the sample to assign the read to.
        """
        # Extract the sequenced barcode from the header of the read
        sequenced_barcode = self.extract_barcode_header(read_pair.forward_read, barcode_length)
        # For training purpose, print the barcode extracted from the read header
        #print("Test: sequenced_barcode: %s\n" % sequenced_barcode)
        # if the barcode exactly matches one of the expected barcodes
        if sequenced_barcode in self.expected.keys():
            # set the assigned sample to the corresponding sample_id
            read_pair.sample = self.expected[sequenced_barcode]
        # if the barcode exactly matches an unexpected barcodes previously resolved to a sample
        elif sequenced_barcode in self.mismatched.keys():
            # set the assigned sample to the corresponding sample_id
            read_pair.sample = self.mismatched[sequenced_barcode]
        # if the barcode exactly matches an unexpected barcodes previously resolved to ambiguous or unmatched ones
        elif sequenced_barcode in self.unmatched:
            # leave the assigned sample to False (unassigned)
            # pass is the Python command that does nothing
            pass
        # if the barcode does not match any barcode observed so far
        else:
            # resolve if it can uniquely be assigned to one expected barcode
            # assigned_barcode will be the sample_id, or False if the barcode is ambiguous or unassigned
            read_pair.sample = self.assign_to_unique_sample(sequenced_barcode)
            # If a unique sample is resolved
            if read_pair.sample:
                # Update the observed barcodes accordingly
                self.mismatched[sequenced_barcode] = read_pair.sample
            else:
                self.unmatched.append(sequenced_barcode)

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
            # check if the sequenced barcode is an approximate match to each of them (1 substitution allowed)
            if ApproxMatch.approx_substitute(sequenced_barcode, expected, 1):
            # The alternative below takes a few more seconds because it initialises additional variable useful for 
            # browsing long sequences to find a smaller sequence. In this case, a direct comparison of an expected
            # barcode and a sequenced barcode of the same length is faster.
            # TODO: for a single mismatch approach, try a Levenshtein approach
            #if len(fuzzysearch.find_near_matches(sequenced_barcode, expected, max_substitutions=1,
            #                                 max_insertions=0, max_deletions=0, max_l_dist=None)):
                matches.append(expected)
        # If the sequenced barcode is an approximate match of a unique expected barcode
        if len(matches) == 1:
            # return the corresponding sample name
            return self.expected[matches[0]]
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
        #print("Test: barcode_header: %s" % self.barcode_pattern.search(read.header_line).group(1))
        return self.barcode_pattern.search(read.header_line).group(1)


class AdapterParser:
    """Parses expected adapter sequence, unique identifier, and read mate to screen for it from a user file.

    This type of objects was given variables and functions to import adapters expected to be present in the reads."""

    def __init__(self, filename):
        """Constructor for AdapterParser"""
        # The name of the file to read the adapters from
        self.filename = filename
        # a named tuple to store the forward and reverse adapter
        self.adapters = None

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
                        'Error: Unexpected value in the "mate" field. Please check your adaptor file (%s) or seek '
                        'advice about this script!' % self.filename)
                    sys.exit(6)
        # Once we finished parsing the adaptor file, we do some sanity checks and if ok, return the adaptor object
        # return an error and exit the script if there are more than 1 adaptor expected in the forward read
        if len(forward_identifiers) > 1:
            print(
                'Error: More than one adaptor sequence is expected in the forward read. Please check your adaptor '
                'file (%s) '
                'or seek advice about this script!' % self.filename)
            sys.exit(7)
            # Return an error
        if len(reverse_identifiers) > 1:
            print(
                'Error: More than one adaptor sequence is expected in the reverse read. Please check your adaptor '
                'file (%s) '
                'or seek advice about this script!' % self.filename)
            sys.exit(8)
        # Return an error if the identifier from forward and reverse read are the same (if we reached this line,
        # we know there is one unique adapter for each mate).
        if forward_identifiers == reverse_identifiers:
            print(
                'Error: Identifiers from forward and reverse are not unique. Please check your adaptor file (%s) or '
                'seek '
                'advice about this script!' % self.filename)
            sys.exit(9)
        # After all the sanity check return the forward and reverse adaptor objects as a tuple
        # Prepare a named tuple structure to cleanly return the two adapters in one variable
        adapter_pair = collections.namedtuple('AdapterPair', ['forward', 'reverse'])
        # Fill in the tuple structure with the adapter objects and return it
        self.adapters = adapter_pair(forward=SeqDataTypes.Adapter(forward_identifiers[0], forward_sequence, 'forward'),
                                     reverse=SeqDataTypes.Adapter(reverse_identifiers[0], reverse_sequence, 'reverse'))


class FastqWriter:
    """Writes fastq files.

    This type of object was given variables and functions to handle multiple file writing streams open simultaneously
    and write the processed reads to the corresponding output fastq file."""

    def __init__(self, basename_raw, samples):
        """Constructor for FastqWriter"""
        # Create a dictionary of file handlers (values) associated to each sample (key)
        self.accepted_outfiles = {}
        self.filtered_outfiles = {}
        self.initialise_fastq_outfiles(basename_raw, samples)

    def initialise_fastq_outfiles(self, basename_raw, samples):
        """Opens file writing streams for each sample to write out the accepted and excluded forward and reverse
        reads, plus two more writing streams for the unassigned forward and reverse reads.

        Note that basename_raw is the name of the input forward read file, not the name of any output file. The original
        basename_raw will be used to name the output files containing the reads unassigned to any sample.

        Args:
            self, samples

        Returns:
            None
        """
        # Prepares a dictionary of pairs of file streams. One pair of file stream corresponds to the forward and
        # reverse mate of a sample, or the unassigned reads. Each pair of file streams will be accessible from the
        # dictionary using the sample_id as the key. The forward file stream will be accessible using "pair.forward".
        outfiles = {}
        # To save pairs of file streams we will use the collections module again
        filestream_pair = collections.namedtuple('FileStreamPair', ['forward', 'reverse'])
        # For each sample name,
        for sample in samples:
            # Open a file stream for the accepted forward read toward a file name sample_pe1.fastq.gz
            self.accepted_outfiles[sample] = filestream_pair(forward=open("%s_pe1.fastq" % sample, 'w'),
                                                    reverse=open("%s_pe2.fastq" % sample, 'w'))
        # For each sample name,
        for sample in samples:
            # Open a file stream for the filtered forward read toward a file name sample_pe1.fastq.gz
            self.filtered_outfiles[sample] = filestream_pair(forward=open("%s_pe1.excluded" % sample, 'w'),
                                                    reverse=open("%s_pe2.excluded" % sample, 'w'))
        # Create an additional pair of file streams for the unassigned reads
        # Note that the sample field of unassigned reads was set to False
        self.filtered_outfiles[False] = filestream_pair(forward=open("%s_1.excluded" % basename_raw, 'w'),
                                               reverse=open("%s_2.excluded" % basename_raw, 'w'))


    def write_accepted_reads(self, read_pair):
        """Write the forward and reverse read in the corresponding filestream_pair of file stream.

            Note that only the forward_read was deconvoluted. The barcode of the reverse read is assumed to be the same.
        Args:
            self, read

        Returns:
            None
        """
        # get the filestream_pair of file stream corresponding to the sample the forward read was deconvoluted to
        filestream_pair = self.accepted_outfiles[read_pair.sample]
        # write the forward read in the forward file
        filestream_pair.forward.write("%s\n" % read_pair.forward_read)
        # write the reverse read in the reverse file
        filestream_pair.reverse.write("%s\n" % read_pair.reverse_read)

    def write_filtered_reads(self, read_pair):
        """Write the forward and reverse read in the corresponding filestream_pair of filestreams.

            Note that only the forward_read was deconvoluted. The barcode of the reverse read is assumed to be the same.
        Args:
            self, read

        Returns:
            None
        """
        # get the filestream_pair of file stream corresponding to the sample the forward read was deconvoluted to
        filestream_pair = self.filtered_outfiles[read_pair.sample]
        # write the forward read in the forward file
        filestream_pair.forward.write("%s\n" % read_pair.forward_read)
        # write the reverse read in the reverse file
        filestream_pair.reverse.write("%s\n" % read_pair.reverse_read)

    def close_files(self):
        """Closes all the writing file streams to save the files.

        Args:
            self

        Returns:
            None
        """
        for pair in self.accepted_outfiles.values():
            pair.forward.close()
            pair.reverse.close()
        for pair in self.filtered_outfiles.values():
            pair.forward.close()
            pair.reverse.close()


def set_basename_from_raw_file(filename):
    """Generates a report filename based on the name of the filename of the raw reads.

    Note that the filename of the raw reads will look different depending on the sequencing centre. For now only the
    BGI case will be implemented to test the AlvMac files.

    Args:
        filename

    Returns:
        None. Sets the self.filename field instead.
    """
    # Extracts the basename of the input file
    basename = os.path.basename(filename)
    #print("Test: basename: %s" % basename)
    # For each know centre, test if the filename matches its typical pattern
    # Test for BGI (Raw_"anything"_1."anything")
    # keep the "Raw_anything" until the character before the _1
    BGI_match = re.search('Raw_.*(?=_1\..*$)', basename)
    test_match = re.search('single.*(?=_pe1.*)|.*thousands(?=_pe1.*)', basename)
    if BGI_match:
        unassigned_basename = BGI_match.group(0)
    # Test for our training files
    # keep the name from the start until the character before the _pe1
    elif test_match:
        unassigned_basename = test_match.group(0)
    # if no pattern matches
    else:
        # return an error and exit
        print(
            "Error: Your data seems to have been generated by an unknown sequencing centre; please seek advice about "
            "this script!")
        sys.exit(2)
    # if we reach here, a basename was defined, complete it with the extension and return the forward and reverse
    # filenames
    return unassigned_basename

