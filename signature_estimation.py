#!/usr/bin/env python

################################################################################
# SET UP
################################################################################
# Load required modules
import sys, os, argparse, pandas as pd, logging
import quadprog, numpy as np
from sklearn.utils import check_array

# Constants
QP = 'QP'
SIG_EST_ALGORITHMS = [ QP ]

# Logging
FORMAT = '%(asctime)s SignatureEstimation %(levelname)-10s: %(message)s'
logging.basicConfig(format=FORMAT)
def get_logger(verbosity=logging.INFO):
    logger = logging.getLogger(__name__)
    logger.setLevel(verbosity)
    return logger

################################################################################
# SIGNATURE ESTIMATION METHODS
################################################################################
# Wrapper
def signature_estimation(M, P, algorithm):
    if algorithm == QP:
        return signature_estimation_qp(M, P)
    else:
        raise NotImplementedError('Algorithm "%s" not implemented.')

# QP
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
# Parse command-line arguments
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mf', '--mutation_counts_file', type=str, required=True)
    parser.add_argument('-a', '--algorithm', type=str, required=False,
                        default=QP, choices=SIG_EST_ALGORITHMS)
    parser.add_argument('-sf', '--signatures_file', type=str, required=True)
    parser.add_argument('-of', '--output_file', type=str, required=True)
    parser.add_argument('-as', '--active_signatures', type=int, required=False,
                        default=[], nargs='*')
    parser.add_argument('-v', '--verbosity', type=str, required=False, default=logging.INFO)
    return parser

# Main
def run( args ):
    # Load the mutation counts
    logger = get_logger(args.verbosity)
    logger.info('[Loading input files]')
    
    mut_df = pd.read_csv(args.mutation_counts_file, sep='\t', index_col=0)
    categories = list(mut_df.columns)
    M = mut_df.values

    logger.info('- Loaded %s x %s mutation count matrix' % M.shape)

    # Load the signatures
    sigs_df = pd.read_csv(args.signatures_file, sep='\t', index_col=0)
    assert(list(sigs_df.columns) == categories)
    sigs = sigs_df.values
    logger.info('- Loaded %s x %s signature matrix' % sigs.shape)
            
    # Restrict to active signatures (if necessary)
    if len(args.active_signatures) > 0:
        logger.info('\t-> Restricting to %s signatures...' % len(args.active_signatures))
        sigs = sigs[np.array(args.active_signatures)-1]
        active_signatures = args.active_signatures
    else:
        active_signatures = sigs_df.index
            
    # Compute the exposures and output to file
    # SignatureEstimation gives the proportion of mutations per signature,
    # so we renormalize by multiplying by the number of mutations per sample
    logger.info('[Computing exposures]')

    n_muts = M.sum(axis=1)
    exp_df = pd.DataFrame(index=mut_df.index, columns=active_signatures,
                          data=signature_estimation(M, sigs, args.algorithm))

    logger.info('- Outputting to file: %s' % args.output_file)
    exp_df.to_csv(args.output_file, sep='\t')

if __name__ == '__main__': run( get_parser().parse_args(sys.argv[1:]) )
