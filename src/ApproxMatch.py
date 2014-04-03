__author__ = 'David Magee, Carolina Correia, and Kevin Rue-Albrecht'
__copyright__ = "Copyright 2014, GPLv2"

"""Empty docstring
"""

# Module arrays allows to handle constrained and efficient array structure
from array import array

def approx_substitute(str1, str2, max_substitutions):
    """Checks that str1 is less than max_substitutions away from str2.

    Note: No insertions or deletions are allowed. Sequences of different length will automatically return FALSE.

    Args:
        str1, str2, max_substitutions

    Returns:
        Boolean. TRUE if str1 is less than ax_substitutions away from str2, FALSE otherwise.
    """
    # Solves a simple scenario which does not require to parse the sequences.
    #if len(str1) != len(str2):
    #    return False
    # If we reach here, we know that the two strings are the same length, therefore len(str1) is synonym to len(str2)
    # For each (index, character) value pair in str1,str2
    for (c1, c2) in zip(str1, str2):
        # remove 1 if the characters from str1 and str2 are different
        max_substitutions -= c1 != c2
        # time saver: if the counter goes below zero, more substitutions than allowed were found
        if max_substitutions < 0:
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
        result = approx_starts_with(small_string, long_string[str2_idx:], max_ins, max_sub)
        if result:
            return result
    return False