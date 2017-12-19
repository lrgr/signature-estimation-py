#!/usr/bin/env python

# Load required modules
import sys, os, pandas as pd, logging, numpy as np, json

################################################################################
# LOGGING
################################################################################
FORMAT = '%(asctime)s %(filename)-15s %(levelname)-10s: %(message)s'
logging.basicConfig(format=FORMAT)

def getLogger(verbosity=logging.INFO):
    logger = logging.getLogger(__name__)
    logger.setLevel(verbosity)
    return logger

################################################################################
# READING FILES
################################################################################
def load_mutation_counts(mutation_counts_file, logger=getLogger()):
    df = pd.read_table(mutation_counts_file, sep='\t', index_col=0)
    samples = list(df.index)
    categories = list(df.columns.values)
    M = df.as_matrix()
    logger.info('- Loaded mutation data in %s samples of %s categories' % M.shape)
    return M, samples, categories

def load_signatures(signature_file, categories, logger=getLogger()):
    # Parse the file
    with open(signature_file, 'r') as IN:
        obj = json.load(IN)
        signatures = obj.get('signatures')
        typeToSignatures = obj.get('typeToSignatures')
        sig_categories = obj.get('categories')
        sig_names = obj.get('signature_names')

    logger.info('- Loaded %s signatures' % len(signatures))

    # Reorder the signatures according to the given categories (if necessary)
    catToIndex = dict(zip(categories, range(len(categories))))
    indices = [ catToIndex[c] for c in sig_categories ]
    signatures = dict( (s, np.array(sig)[indices]) for s, sig in signatures.items() )
    P = np.array([ signatures[s] for s in sig_names ])

    return P, sig_names, typeToSignatures
