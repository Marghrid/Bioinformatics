# Viterbi Algorithm

viterbi.py implements the algorithms asked in section III of the problem statement.

The package numpy is required to run the script.

Usage:
python3 viterbi.py [--prob] [--post=n] S

where S is the sequence we want to analyse.

Using --prob calculates the probability of observing the sequence (P[S])
Using --post=n calculates the posterior probabilities for state n (P[pi_n|S])

Example:
python3 viterbi.py --prob CATGCGGGTTATAAC

## Implementation Details

We use matrix operations (through the use of the numpy package) to facilitate the implementation of the
Viterbi algorithm.

We also normalize (2-norm) the vector m_t (which would correspond to the probabilities of each sequence up to time t)
at the end of each timestep in order to avoid problems with underflowing. This is correct because the 2-norm
normalization is monotonous.
