#!/usr/bin/env python

# Load required modules
import sys, os, argparse, json, logging, numpy as np

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

# Load the signatures file
with open(args.input_file, 'r') as IN:
    arrs = [ l.rstrip().split('\t') for l in IN ]
    header = arrs.pop(0)
    signature_names = header[3:]
    categories = [ arr[2] for arr in arrs ]
    signature_matrix = np.array([ [ float(p) for p in arr[3:] ] for arr in arrs ]).T
    signatures = dict(zip(signature_names, signature_matrix.tolist()))

logger.info('Loaded %s x %s signature matrix...' % signature_matrix.shape)

# Output to file
with open(args.output_file, 'w') as OUT:
    output = dict(params=vars(args), signature_names=signature_names,
                  signatures=signatures, categories=categories)
    json.dump(output, OUT)
