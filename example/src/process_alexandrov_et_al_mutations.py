#!/usr/bin/env python

# Load required modules
import sys, os, argparse, pandas as pd, logging, numpy as np

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', type=str, required=True)
parser.add_argument('-o', '--output_file', type=str, required=True)
parser.add_argument('-v', '--verbosity', type=int, default=logging.INFO, required=False)
args = parser.parse_args(sys.argv[1:])

# Set up logger
logger = logging.getLogger(__name__)
logger.setLevel(args.verbosity)

# Helpers for parsing categories into substitution, left flanking,
# and right flanking
def sub(c): return c.split('[')[1].split(']')[0]
def lf(c): return c.split('[')[0]
def rf(c): return c.split(']')[-1]

# Load the input file and convert to Pandas DataFrame
with open(args.input_file, 'r') as IN:
    arrs = [ l.rstrip('\n').split('\t') for l in IN ]
    header = arrs.pop(0)
    samples = header[1:]
    cats = [ arr[0] for arr in arrs ]
    counts = np.array([ [int(c) for c in arr[1:] ] for arr in arrs ]).T

    # Sort the categories in the standard way
    sorted_cats = sorted(cats, key=lambda c: (sub(c), lf(c), rf(c)))
    cat_idx = [ cats.index(c) for c in sorted_cats ]
    counts = counts[:, cat_idx]

df = pd.DataFrame(data=counts, columns=sorted_cats, index=samples)

logger.info('Loaded dataset with...')
logger.info('- %s samples' % len(samples))
logger.info('- %s categories' % len(cats))
logger.info('- %s mutations' % counts.sum())

# Output to file
df.to_csv(args.output_file, sep='\t', index=True)
