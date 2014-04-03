__author__ = 'David Magee, Carolina Correia, and Kevin Rue-Albrecht'
__copyright__ = "Copyright 2014, GPLv2"

"""Empty docstring
"""

# Module fuzzysearch allows approximate matching of a substring within a larger string
import fuzzysearch
# Custom Module which contains a few functions for approximate matching
import ApproxMatch

class Read:
    """Sequenced read and associated data"""

    def __init__(self, header_line, sequence_line, separator_line, quality_line):
        """Constructor for Read"""
        self.header_line = header_line
        self.sequence_line = sequence_line
        self.separator_line = separator_line
        self.quality_line = quality_line
        # If acceptable quality then set to True, otherwise set to False
        # Defaults to True (if no quality check is done, quality is assumed to be good)
        self.quality_status = True
        # If adapter is present then set to True, otherwise set to False
        # Defaults to False (if no adapter check is done, it is assumed to be absent)
        self.adapter_present = False

    def trim(self, start, stop):
        """Crops the read sequence and associated quality_line keeping positions start to stop.
        Be careful to give the positions in Python position numbering (starts at 0). The pain of calculating the proper
        coordinates before calling the function will save a lot of time avoiding re-calculating them every time the
        function is called.

        Args:
            self, start, stop

        Returns:
            None
        """
        # Trims the sequence line to the desired window
        self.sequence_line = self.sequence_line[start:stop]
        # Trims the quality line to the same window
        self.quality_line = self.quality_line[start:stop]

    def define_quality_status(self, threshold, max_bases, quality_alphabet, length_alphabet, ascii_counts):
        """Sets the quality_status attribute according to whether stricly more than percentage of bases in the read are
         strictly below phred.
         Be careful to give a threshold that is the desired minimal Phred score -1. The reason can be demonstrated with
         the extreme example where the read has the minimal number of bases at the minimal score or higher. The pain of
         calculating the proper coordinates before calling the function will save a lot of time avoiding re-calculating
         them every time the function is called.

        :rtype : None
        Args:
            self, phred, percentage

        Returns:
            None
        """
        # Count how many of each symbol are present
        # this will crash for unexpected characters in the quality line
        # the range was defined/hard-coded in the main script (Processing_paired_reads.py) and may need to be updated 
        # with the technology. The expected range is limited to avoid many ASCII with 0 ascii_counts.
        for symbol in self.quality_line:
            ascii_counts[symbol] += 1
        # Search the smallest ascii with non-zero count
        ascii_index = 0
        while ascii_index < length_alphabet and ascii_counts[quality_alphabet[ascii_index]]==0 :
            ascii_index += 1
        # Return an error and exit if no count was found
        if ascii_index == len(ascii_counts):
            print("Error: No count found for any ASCII symbol in quality line: ", ascii_counts)
            sys.exit(3)
        #print("Test: ascii_counts: ", ascii_counts)
        #print("Test: ascii_counts(sum of values): ", sum(ascii_counts.values()))
        # Given max_bases allowed under Phred threshold, add up the poor quality bases until max_bases is reached
        # and if the corresponding Phred (percentile) is below the threshold, then it means max_bases is exceeded
        while max_bases > 0:
            #print("Test: max_bases: ", max_bases)
            #print("Test: symbol: \"%s\" ascii_counts: %i" % (quality_alphabet[ascii_index], ascii_counts[quality_alphabet[ascii_index]]))
            if ascii_counts[quality_alphabet[ascii_index]] > max_bases:
                break
            max_bases -= ascii_counts[quality_alphabet[ascii_index]]
            ascii_index += 1
            while ascii_index < length_alphabet and ascii_counts[quality_alphabet[ascii_index]]==0 :
                ascii_index += 1
        #print("Test:  max_bases_left = %i, Phred_score = %i, Ascii_value=%i\n" % (max_bases, ascii_index, ascii_index+ord(quality_alphabet[0])))
        if ascii_index < threshold:
            self.quality_status = False
        return

    def define_adapter_presence_substitutions_only(self, adapter, max_substitutions):
        """Sets the adapter_absent attribute according to whether a match is found with a number of substitutions
        less or equal to
        edit_threshold.

        :rtype : None
        Args:
            self, adapter, edit_threshold

        Returns:
            None

        Args:
            self, adapter, max_substitutions

        Returns:
            None
        """
        # Solves the straightforward case where the adapter is exactly present in the read
        #if adapter in self.sequence_line:
        #    self.adapter_present = True
        #    return
        # Otherwise, look for an approximate match of the adapter, less than max_substitutions different
        # First do a preliminary (hopefully faster) check whether the adapter is present as an approximate match in
        # the read sequence within a Levenshtein distance equal to the allowed number of substitutions
        # Two sequences less than N substitutions apart are automatically less than a N edits apart (the reciprocal
        # is not true). If the Levenshtein code is faster, then it will save a lot of time pre-filtering adapters
        # within the Levenshtein distance, which are much more likely to contain adapters. The slower code checking
        # the mismatch distance will then be called on these pre-filtered reads to confirm whether it is an actual
        # substituted match or if the Levensthein match involved insertions and deletions
        # If 1 or more approximate matches of the adapter were found within a Levenshtein distance of max_substitutions
        if fuzzysearch.susbstitutions_only.has_near_match_substitutions_ngrams(adapter, self.sequence_line, max_substitutions):
            self.adapter_present = True
            return
        # The simple fact of arriving here proves that no match was found, therefore leave the adapter
        # presence to False and simply leave the function

    def __str__(self):
        """This function defines the string representation displayed when calling the print function on a Read object"""
        # Uses the new-line separator of the operating system running this script to separate the different elements of
        # the read on different lines at the screen.
        return "\n".join([self.header_line, self.sequence_line, self.separator_line, self.quality_line])
        # returns the header line, sequence line, separator line and quality line of the read on separate lines


class ReadPair:
    """Sequenced read mates and associated data."""

    def __init__(self, forward_read, reverse_read):
        """Constructor for ReadPair"""
        self.forward_read = forward_read
        self.reverse_read = reverse_read
        # If the read was assigned to a sample then set to the sample_id, otherwise set to False
        # Default to False (unassigned read)
        self.sample = False

    def trim(self, start, stop):
        """Crops the read sequence and associated quality_line keeping positions start to stop.
        Be careful to give the positions in Python position numbering (starts at 0). The pain of calculating the proper
        coordinates before calling the function will save a lot of time avoiding re-calculating them every time the
        function is called.

        Args:
            self, start, stop

        Returns:
            None
        """
        self.forward_read.trim(start, stop)
        self.reverse_read.trim(start, stop)


    def __str__(self):
        """String representation of the read pair.

        Args:


        Returns:
            None
        """
        return "%s\n%s" % (self.forward_read, self.reverse_read)


class Adapter:
    """Adapter sequence and associated data"""

    def __init__(self, identifier, sequence, mate):
        """Constructor for Adapter"""
        self.identifier = identifier
        self.sequence = sequence
        self.mate = mate

    def __str__(self):
        """This function defines a string representation displayed when calling the print function on an Adapter
        object"""
        # Uses the tab separator to separate the different elements of the Adapter
        return '\t'.join([self.sequence, self.identifier, self.mate])
        # returns the header line, sequence line, separator line and quality line of the read on separate lines

