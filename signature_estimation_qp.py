#!/usr/bin/env python

# Load required modules
import sys, os, logging, pandas as pd
from i_o import load_mutation_counts, load_signatures, getLogger

################################################################################
# SIGNATURE ESTIMATION METHODS
################################################################################
import quadprog, numpy as np
from sklearn.utils import check_array
def signature_estimation_qp(M, P):
    """
    Estimate exposures with quadratic programming.

    Inputs:
    - M: (N x L) mutation count matrix (N=# of samples, L=# of categories)
    - P: (K x L) mutation signature matrix (K=# of signatures)

    Outputs:
    - E: (N x K) exposure matrix

    Much of this is taken/adapted from the SignatureEstimation R package:
    https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/software/signatureestimation/SignatureEstimation.pdf
    """
    # Do some checks
    P = check_array(P)
    M = check_array(M)

    # Normalize M and transpose to match the SignatureEsimation package
    M = M/M.sum(axis=1)[:, None]
    M = M.T
    P = P.T

    # K: number of signatures
    K = P.shape[1]
    # N: number of samples
    N = M.shape[1]
    # G: matrix appearing in the quatric programming objective function
    G = P.T.dot(P)
    # C: matrix constraints under which we want to minimize the quatric programming objective function.
    C = np.hstack((np.ones((K, 1), dtype=np.float64), np.eye(K, dtype=np.float64)))
    # b: vector containing the values of b_0.
    b = np.array([1.] + [0.] * K, dtype=np.float64)
    # d: vector appearing in the quadratic programming objective function as a^T
    D = M.T.dot(P)

    # Solve quadratic programming problems
    exposures = np.array([ quadprog.solve_qp(G, d, C, b, meq=1)[0] for d in D ])

    # Some exposure values may be negative due to numerical issues,
    # but very close to 0. Change these neagtive values to zero and renormalize.
    exposures[exposures < 0] = 0
    exposures = exposures/exposures.sum(axis=1)[:, None]

    return exposures

################################################################################
# MAIN
################################################################################
def get_parser():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-mf', '--mutation_counts_file', type=str, required=True)
    parser.add_argument('-sf', '--signature_file', type=str, required=True)
    parser.add_argument('-s', '--signatures', type=str, required=False, default=[], nargs='*')
    parser.add_argument('-of', '--output_file', type=str, required=True)
    parser.add_argument('-v', '--verbosity', type=int, default=logging.INFO, required=False)
    return parser

def run(args):
    # Setup logger
    logger = getLogger(args.verbosity)

    # Load the mutation counts and signatures
    logger.info('* Loading mutation data and signatures')
    M, samples, categories = load_mutation_counts(args.mutation_counts_file, logger)
    P, sigs, typeToSignatures = load_signatures(args.signature_file, categories, logger)

    # Restrict to given signatures
    if len(args.signatures) > 0:
        logger.info('* Restricting to %s signatures...' % len(args.signatures))
        sig_indices = [sigs.index(s) for s in args.signatures ]
        P = P[sig_indices]
        sigs = args.signatures

    # Run QP and output to file
    E = signature_estimation_qp(M, P)
    df = pd.DataFrame(E, index=samples, columns=sigs)
    df.to_csv(args.output_file, sep='\t')

if __name__ == '__main__':
    run( get_parser().parse_args(sys.argv[1:]))
