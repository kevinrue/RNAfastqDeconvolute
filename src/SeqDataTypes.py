__author__ = 'David Magee, Carolina Correia, and Kevin Rue-Albrecht'
__copyright__ = "Copyright 2014, GPLv2"

"""Empty docstring
"""

# Module os allows to access OS-specific line separator.
import os
# Module scipy allows to calculate quantiles, useful to define the PHRED score below which a given percentage of the
# bases in a read are found.
import scipy
# Module fuzzysearch allows to find approximate matches, useful to detect adaptor sequences in reads with possible
# sequencing errors.
import fuzzysearch
# Module array allows to efficiently store constrained values, useful for sequence matching.
from array import array


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
        self.adapter_present = True

    def trim(self, left, right):
        """Crops the read sequence and associated quality_line keeping positions five_start to three_end.
        Be careful to give the positions in Python position numbering (starts at 0). The pain of calculating the proper
        coordinates before calling the function will save a lot of time avoiding re-calculating them every time the
        function is called.

        :rtype : None
        Args:
            self, left, right

        Returns:
            None
        """
        # Trims the sequence line to the desired window
        self.sequence_line = self.sequence_line[left:right]
        # Trims the quality line to the same window
        self.quality_line = self.quality_line[left:right]

    def define_quality_status(self, threshold, percentage):
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
        # Converts all the ascii caracters in the quality line into their ascii code
        quality_ascii = [ord(c) for c in self.quality_line]
        # Defines whether the read is acceptable (True) or not (False)
        self.quality_status = scipy.percentile(quality_ascii, percentage) > threshold
        # For training purpose: prints the value of the percentile used for threshold TODO remove in final code
        print('percentile ({0}%): {1:.3f}'.format(percentage, scipy.percentile(quality_ascii, percentage)))


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
        # for each subsequence of the sequence_line of length identical to the adapter (last starting position is
        # length of the adapter before the end of the read)
        for index in range(len(self.sequence_line) - len(adapter) + 1):
            # check if the adapter is less than max_substitutions away from the subsequence
            result = approx_substitute(adapter, self.sequence_line[index:index + len(adapter)], max_substitutions)
            # If a match was found (result is TRUE)
            if result:
                # return TRUE (and stop the function here)
                return result
        # If no substring of sequence_line returned TRUE, the function will eventually leave the loop above
        # The simple fact of arriving here proves that no match was found, therefore return FALSE
        return False

    def define_adapter_presence_levenshtein(self, adapter, edit_threshold):
        """Set the adapter_absent attribute according to whether a match is found with an edit distance less or equal to
        edit_threshold.

        :rtype : None
        Args:
            self, adapter, edit_threshold

        Returns:
            None
        """
        self.adapter_present = len(
            fuzzysearch.find_near_matches_with_ngrams(adapter, self.sequence_line, edit_threshold)) > 0

    def __str__(self):
        """This function defines the string representation displayed when calling the print function on a Read object"""
        # Uses the new-line separator of the operating system running this script to separate the different elements of
        # the read on different lines at the screen.
        return os.linesep.join([self.header_line, self.sequence_line, self.separator_line, self.quality_line])
        # returns the header line, sequence line, separator line and quality line of the read on separate lines


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


def approx_substitute(str1, str2, max_substitutions):
    """Checks that str1 is less than max_substitutions away from str2.

    Note: No insertions or deletions are allowed. Sequences of different length will automatically return FALSE.

    Args:
        str1, str2, max_substitutions

    Returns:
        Boolean. TRUE if str1 is less than ax_substitutions away from str2, FALSE otherwise.
    """
    # Solves a simple scenario which does not require to parse the sequences.
    if len(str1) != len(str2):
        return False
    # If we reach here, we know that the two strings are the same length, therefore len(str1) is synonym to len(str2)
    # Initialise a counter of substitutions between the two strings
    substitutions = 0
    # For each (index, character) value pair in str1
    for (index, str1_char) in enumerate(str1):
        # add 1 if the characters from str1 and str2 are different, add 0 otherwise
        substitutions += (0 if str1_char == str2[index] else 1)
        # time saver: if the counter exceeds max_substitutions at some stage, don't bother checking the rest...
        if substitutions > max_substitutions:
            # ... just return FALSE
            return False
    # If max_substitutions is never reached, the function will eventually leave the loop above
    # The simple fact of arriving here proves that str1 is less than max_substitutions away from str2,
    # therefore return TRUE
    return True


def approx_starts_with(small_string, long_string, max_substitutions, max_insertions):
    """Checks if small_string can be found at the start of long_string, given the restrictions on substitutions and
    insertion allowances.

    Args:
        small_string, long_string, max_substitutions, max_insertions

    Returns:
        Boolean. TRUE if small_string is found in long_string within restrictions, FALSE otherwise.
    """
    # check a simple case which is obviously impossible
    if not len(small_string) <= len(long_string):
        return False
    # initialises to 0 a list of scores for matches with 0, 1, .., max_insertions mismatches
    scores = array('I', [0] * (max_insertions + 1))
    # initialises to 0 a temporary version of the scores above
    new_scores = scores[:]
    # For each (index, character) value pair in small_string
    for (str1_idx, char1) in enumerate(small_string):
        # Excessive value to force the first iteration to use the other value  calculated in the inner loop
        prev_score = len(long_string)
        # For (index, character) value pairs in long_string, from the current position in small_string (
        # str1_index) to that
        # position plus the maximal number of insertions allowed (time saver to prevent larger parsing).
        for (n_insertions, char2) in enumerate(
                long_string[str1_idx:max_insertions + str1_idx + 1]
        ):
            # updates the scores above, adds 1 to the situations with a substitution, and leave the score as is
            # if equal. The lowest score will describe the closest match of small_string in long_string:
            # the value of the score will be the number of substitutions
            # the index of the score will be the number of insertions
            new_scores[n_insertions] = prev_score = min(scores[n_insertions] + (0 if char1 == char2 else 1),
                                                        prev_score)

        # swap scores <-> new_scores, to save the current scores in "scores" and re-initialise "new_scores" for the
        # next iteration of the loop
        scores, new_scores = new_scores, scores
    # return if the lowest score (= number of substitutions) is below the allowed threshold
    # the check for insertions below the allowed threshold is implied as scores are only calculated for situations
    # between 0 and the maximal number of insertions anyway.
    return min(scores) <= max_substitutions


def approx_within(small_string, long_string, max_ins, max_sub):
    """Check if it is possible to find small_string within long_string given limitations

    Args:
        target, reference, max_substitutions, max_insertions

    Returns:
        Boolean. TRUE if target is found in reference and within restrictions, FALSE otherwise.
    """
    for str2_idx in range(len(long_string) - len(small_string) + 1):
        print(str2_idx)
        result = approxStartsWith(small_string, long_string[str2_idx:], max_ins, max_sub)
        if result:
            return result
    return False