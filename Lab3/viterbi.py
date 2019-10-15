#!/bin/python

import argparse
import numpy as np

hidden_states = ["A", "T", "G", "C"]
states = ["1", "2", "3"]

starting_distribution = np.ones(len(states)) / len(states)

transitions_probabilities = np.array([[.6, .4, 0],
                                      [.25, .5, .25],
                                      [.25, .25, .5]])

emission_probabilities = np.array([[.4, .3, .3, 0],
                                   [.1, .1, .4, .4],
                                   [.4, .3, 0, .3]])


def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Computes the most likely sequence of states, given a sequence of DNA.')
    parser.add_argument('sequence', metavar='S', type=str, help='the sequence')
    parser.add_argument('--prob', action='store_true', help='compute the probability of the sequence of states')
    args = parser.parse_args()
    S = args.sequence

    print("S    =", S)

    traces = np.zeros((len(S), len(states)), dtype=int)
    pi = np.zeros((len(S),), dtype=int)

    # Forward pass
    m = np.diag(emission_probabilities[:, hidden_states.index(S[0])]).dot(starting_distribution)
    for t in range(1, len(S)):
        m_new = normalize(np.diag(emission_probabilities[:, hidden_states.index(S[t])]).dot(
            np.max(transitions_probabilities.T.dot(np.diag(m)), axis=1)))
        traces[t] = np.argmax(transitions_probabilities.T.dot(np.diag(m)), axis=1)
        m = m_new

    # Find max at the end
    pi[len(S) - 1] = np.argmax(m)

    # Traceback
    for t in range(len(S) - 2, -1, -1):
        pi[t] = traces[t + 1][pi[t + 1]]

    # Output optimal path
    print('π*   = ', end='')
    for s in pi:
        print(states[s], end='')
    print()

    # Compute P[S] using the forward algorithm
    if args.prob:
        a = np.diag(emission_probabilities[:, hidden_states.index(S[0])]).dot(starting_distribution)
        for t in range(1, len(S)):
            a = np.diag(emission_probabilities[:, hidden_states.index(S[t])]).dot(transitions_probabilities.T).dot(a)

        print("P[S] =", np.sum(a))
