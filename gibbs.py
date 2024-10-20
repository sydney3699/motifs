"""gibbs sampling motif finder"""
import numpy as np

"""
GibbsMotifFinder(DNA, k-length)
    1. Randomly pick k-length subsequences from each DNA sequence as initial Motifs
    2. Repeat for a set number of iterations (e.g., 10,000) or until convergence:
        3. Randomly select one of the DNA sequences to leave out (DNA_i)
        4. Construct a PWM from all Motifs except the one from the excluded sequence (Motif_i)
        5. Use the PWM to score all k-length subsequences in the excluded sequence (DNA_i)
        6. Select a new position m for the motif in DNA_i **probabilistically**, based on the PWM scores
        7. Update the motif for DNA_i with the new k-mer at position m
    8. Return the final PFM (Position Frequency Matrix)

also looking for promoter sequences? has CDS designation and TATA box, forward or reverse strand
"""

def score_kmer():
      """score kmer, take out of log2 space with np.exp2, then add to scores list"""


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


    # Iteration:
        # 1. choose a sequence at random (from initial seq list)
            
        # 2. recalculate PWM, omitting the chosen seq (seq_i)
            
        # 3. score each possible motif in chosen seq against PWM
            # each score is sum of PWM values at each position(i) & nucleotide(j) in motif
            # score = sum(w(i,j)) for every posible motif in seq_i


        # 4. choose new motif based on score distribution (graph is position v. distribution?)
            

