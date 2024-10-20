"""gibbs sampling motif finder"""
import numpy as np
from motif_ops import build_pfm, build_pwm
from seq_ops import get_seq, reverse_complement

"""
GibbsMotifFinder(DNA, k-length)
    1. Randomly pick k-length subsequences from each DNA sequence as initial Motifs

    2. Repeat for a set number of iterations (e.g., 10,000) or until convergence:

        3. Randomly select one of the DNA sequences to leave out (DNA_i)
           
        4. Construct a PWM from all Motifs except the one from the excluded sequence (Motif_i)
            - Build a PWM using the current motifs excluding the motif from the randomly selected sequence.
            - This allows for computing the expected pattern without bias from that sequence.

        5. Use the PWM to score all k-length subsequences in the excluded sequence (DNA_i)
            - For each k-mer in DNA_i, calculate a score using the PWM.
            - Convert these scores into a probability distribution (higher-scoring k-mers have higher probability).

        6. Select a new position m for the motif in DNA_i probabilistically, based on the PWM scores
            - Choose a new k-mer from DNA_i based on the score distribution. 

        7. Update the motif for DNA_i with the new k-mer at position m
            
    8. Return the final PFM (Position Frequency Matrix)

looking for promoter sequence; has CDS feature type and "AGGAGG" (S-D seq) present, forward or reverse strand
    - if S-D in promoter and promoter has CDS, can be added to list of promoters (to be used for gibbs?)

"""


""" need to figure out where/how to use score_kmer, score_sequence, get_seq, reverse_complement """


def GibbsMotifFinder(DNA, k, max_iter=10000):
    N = len(DNA)  # Number of DNA sequences
    motifs = [random_kmer(seq, k) for seq in DNA]  # Random initial motifs

    for _ in range(max_iter):
        i = random.randint(0, N - 1)  # Randomly select one sequence
        excluded_motifs = [motif for j, motif in enumerate(motifs) if j != i]
        
        # Build PWM from motifs except the one in sequence i
        pwm = build_pwm(excluded_motifs, k)
        
        # Select a new motif probabilistically for sequence i
        new_motif = select_kmer_probabilistically(DNA[i], pwm, k)
        
        # Update the motif for sequence i
        motifs[i] = new_motif

    # Convert final motifs into a Position Frequency Matrix (PFM)
    return build_pfm(motifs, k)

def random_kmer(sequence, k):
    """Pick a random k-mer from the sequence."""
    start = random.randint(0, len(sequence) - k)
    return sequence[start:start + k]

def select_kmer_probabilistically(sequence, pwm, k):
    """Select a k-mer probabilistically based on PWM scores."""
    scores = []
    nuc_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        score = 1
        for pos, nucleotide in enumerate(kmer):
            score *= pwm[nuc_index[nucleotide], pos]  # Multiply probabilities
        scores.append(score)

    # Normalize scores to get probabilities
    probs = [score / sum(scores) for score in scores]
    
    # Select a k-mer based on the probability distribution
    return random.choices([sequence[i:i + k] for i in range(len(sequence) - k + 1)], probs)[0]

def score_kmer(seq, pwm):
    """Function to score a kmer with a pwm
        kmer length is expected to be the same as pwm length

    Args:
        seq(str): kmer to score
        pwm (numpy array): pwm for scoring

    Yields:
        score (float): PWM score for kmer
    """
    
    # Initialize score to 0
    score = 0
    
    if len(seq) != len(pwm[0]):
        raise ValueError('K-mer and PWM are different lengths!')
    
    # Translator for DNA to numeric indices
    dna_to_index = str.maketrans('ACGT', '0123')
   
    # Iterate across kmer and sum log likelihoods
    for j, val in enumerate(list(seq.translate(dna_to_index)), 0):
        score += pwm[int(val), j]

    # Return score
    return score

"""
def  GibbsMotifFinder(seqs, k, seed=42):
    # identify a common motif of length k that occurs once in each seq
    # use Gibbs sampling to iteratively improve guesses for the location of this motif across all sequences.

    # Initialization:
        # 1.randomly choose k-length subsequences from each DNA sequence as initial 'motifs'
        motifs = {}
        for seq in seqs:
            start = np.random.randint(0, len(seq) - k + 1)
            motif = seq[start:start + k]
            seq_without_motif = seq[:start] + seq[start + k:]
            motifs[seq_without_motif] = motif
            
        # 2. build PFM:
            # calculate the freq of every nucleotide at each positon within the motif
            # q(i,j) = frequency of nucleotide i in motif position j      
        nuc_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        pfm = np.zeros((4, k))
        for motif in motifs:
            for position, letter in enumerate(motif):
                row = nuc_index[letter]
                pfm[row, position] += 1
        pfm /= len(motifs)

        # 3. calculate the freq of of every nucleotide in non-motif positions as background
            # bg(j) = frequency of nucleotide j in background
            # (can use 0.25 for each, then pseudocount p = 0.25)
        max_len = max(len(seq) for seq in seq_without_motif)
        bg_pfm = np.zeros((4, max_len))
        total_nucs = 0
        for seq in seq_without_motif:
            for position, letter in enumerate(seq):
                row = nuc_index[letter]
                bg_pfm[row, position] += 1
                total_nucs += 1
        nuc_sums = np.sum(bg_pfm, axis=1)
        bg_freqs = nuc_sums / total_nucs
        bg_pfm /= len(seq_without_motif)
         

        # 4. convert to/build PWM: 
            # w(i,j) = PWM(i,j) formula with either p or bg(j) values
        
        # pwm calculation with bg 


        # pwm calculation with p  
        # Initialize an empty numpy array
        pwm = np.zeros((4,len(pfm[0])), dtype=float)
    
        # Calculate the sums of each column
        sums = np.sum(pfm, axis = 0)
        p = .25
        bg = .25
    
        # For each position in the matrix, apply the PWM formula (using pseudocount)
        for i in range(4):
                for j in range(len(pwm[0])):
                        pwm[i,j] = np.log2((pfm[i,j] + p) / (sums[j]+p*4)) - np.log2(bg)
"""

    # Iteration:
        # 1. choose a sequence at random (from initial seq list)
            
        # 2. recalculate PWM, omitting the chosen seq (seq_i)
            
        # 3. score each possible motif in chosen seq against PWM
            # each score is sum of PWM values at each position(i) & nucleotide(j) in motif
            # score = sum(w(i,j)) for every posible motif in seq_i


        # 4. choose new motif based on score distribution (graph is position v. distribution?)
            

