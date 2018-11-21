#!/usr/bin/env python

################################################################################
# SET UP
################################################################################
# Load required modules
import sys, os, argparse, pandas as pd, logging
import quadprog, numpy as np
from sklearn.utils import check_array
from signature_estimation import signature_estimation

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
    parser.add_argument('-s', '--sample', type=str, required=True)
    parser.add_argument('-ni', '--num_iterations', type=int, required=True)
    parser.add_argument('-v', '--verbosity', type=str, required=False, default=logging.INFO)
    return parser

# Main
def run( args ):
    # Load the mutation counts
    logger = get_logger(args.verbosity)
    logger.info('[Loading input files]')

    mut_df = pd.read_csv(args.mutation_counts_file, sep='\t', index_col=0)
    categories = list(mut_df.columns)
    # M = mut_df.values

    # logger.info('- Loaded %s x %s mutation count matrix' % M.shape)

    # Load the signatures
    sigs_df = pd.read_csv(args.signatures_file, sep='\t', index_col=0)[categories]
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
    M = mut_df.loc[args.sample]
    n_muts = int(M.sum())
    print(n_muts)
    output = []
    print(M)
    print(sigs)
    for i in range(0, args.num_iterations):
        # sample with replacement
        print(M.shape)
        bootstrap_M = M.iloc[np.random.randint(0, M.shape[0], size=n_muts)].index.value_counts()
        bootstrap_M = bootstrap_M.reindex(index=M.index, fill_value=0)
        # estimation exposures
        exposures = signature_estimation(bootstrap_M.values.reshape(1,-1), sigs, args.algorithm)
        # append to output
        output.append(pd.DataFrame(exposures*n_muts))
    # concat index=range(0, args.num_iterations)
    df = pd.concat(output)
    df.index = range(0, args.num_iterations)
    df.columns = active_signatures
    print(df)

    logger.info('- Outputting to file: %s' % args.output_file)
    df.to_csv(args.output_file, sep='\t')

if __name__ == '__main__': run( get_parser().parse_args(sys.argv[1:]) )
