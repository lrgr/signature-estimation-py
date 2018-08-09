#!/usr/bin/env python

# Load required modules
import matplotlib
matplotlib.use('Agg')
import sys, os, argparse, pandas as pd
import matplotlib.pyplot as plt, seaborn as sns

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', type=str, required=True, help='Exposures file')
parser.add_argument('-t', '--title', type=str, required=True)
parser.add_argument('-o', '--output_file', type=str, required=True)
args = parser.parse_args(sys.argv[1:])

# Load the input exposures and expand for our purposes
exp_df = pd.read_csv(args.input_file, sep='\t', index_col=0)
items = []
sig_names = [ s.split(' ')[-1] if type(s) == type('') else s for s in exp_df.columns ]
for patient, r in exp_df.iterrows():
    n_mutations = sum( r[sig] for sig in exp_df.columns )
    for i, sig in enumerate(exp_df.columns):
        items.append({
            "Patient": patient,
            "Signature": sig_names[i],
            "Exposure": r[sig],
            "Prevalence": r[sig]/n_mutations
        })

df = pd.DataFrame(items)
        
# Create plots and output to file
sns.set_style('whitegrid')
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14,5))
sns.boxplot(x="Signature", y="Exposure", order=sig_names, data=df, ax=ax1)
sns.boxplot(x="Signature", y="Prevalence", order=sig_names, data=df, ax=ax2)
plt.suptitle(args.title, y=1)

# Save to file
plt.tight_layout()
plt.savefig(args.output_file)
