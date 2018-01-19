#!/usr/bin/env python

# Load required modules
import sys, os, argparse, pandas as pd, logging, numpy as np

# Load our modules
this_dir = os.path.dirname(__file__)
sys.path.append(os.path.normpath(os.path.join(this_dir, '../')))
from i_o import getLogger

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', type=str, required=True)
parser.add_argument('-o', '--output_file', type=str, required=True)
parser.add_argument('-v', '--verbosity', type=int, default=logging.INFO, required=False)
args = parser.parse_args(sys.argv[1:])

# Set up logger
logger = getLogger(args.verbosity)

# Load the input file and convert to Pandas DataFrame
with open(args.input_file, 'r') as IN:
    arrs = [ l.rstrip('\n').split('\t') for l in IN ]
    header = arrs.pop(0)
    samples = header[1:]
    categories = [ arr[0] for arr in arrs ]
    counts = np.array([ [int(c) for c in arr[1:] ] for arr in arrs ]).T

df = pd.DataFrame(data=counts, columns=categories, index=samples)

logger.info('Loaded dataset with...')
logger.info('- %s samples' % len(samples))
logger.info('- %s categories' % len(categories))
logger.info('- %s mutations' % counts.sum())

# Output to file
df.to_csv(args.output_file, sep='\t', index=True)
