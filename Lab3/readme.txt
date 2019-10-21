viterbi.py implements the algorithms asked in section III of the problem statement.

The package numpy is required to run the script.

Usage:
python3 viterbi.py [--prob] S

where S is the sequence we want to analyse.

Using --prob calculates the probability of observing the sequence (P[S])

Example:
python3 viterbi.py --prob CATGCGGGTTATAAC