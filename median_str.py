"""median string problem"""

def hamming_distance(s1, s2):
    """Calculate the Hamming distance between two strings."""
    
    s1 = s1.upper() # make sure all letters are uppercase
    s2 = s2.upper()
    minlength = min(len(s1), len(s2)) # use shortest len just in case
    count = 0
    for i in range(minlength): # iterate through the min length for both strings
        if s1[i] != s2[i]:
            count += 1 # add 1 if different
    return count
    

def total_distance(pattern, dna_list):
    """Calculate the total distance between a pattern and a list of DNA sequences."""
    # TODO: Implement the total distance calculation

    # 
    pass


def branch_and_bound(prefix, k, dna_list, best_distance):
    """
    Recursive branch and bound algorithm for median string search.
    
    bfs uses queue (first in first out), dfs uses stack (last in last out), both are O(n)

    Args:
        prefix: Current prefix of the pattern being built.
        k: Remaining length to complete the pattern.
        dna_list: List of DNA sequences.
        best_distance: Current best distance found.
    
    Returns:
        A tuple containing the best pattern found and its distance.
    """
    # TODO: Implement the branch and bound algorithm
    pass


def median_string(dna_list, k):
    """
    Find the median string of length k for a list of DNA sequences using branch and bound.
    
    Args:
        dna_list: List of DNA sequences.
        k: Length of the median string to find.
    
    Returns:
        A tuple containing the median string of length k and its total distance.
    """
    # TODO: Implement the median string search using the branch_and_bound function
    pass


def validate_input(dna_list, k):
    """Validate the input DNA sequences and k value."""
    # TODO: Implement input validation
    # Raise appropriate exceptions for invalid inputs
    pass