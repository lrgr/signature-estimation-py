# Signature Estimation in Python
<img src='https://travis-ci.com/lrgr/signature-estimation-py.svg?token=xpopk4qvQzXty9qXHH3S&branch=master'>

Adapting the [Signature Estimation R](https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/index.cgi#signatureestimation) package for Python. Currently, this package only implements the quadratic programming (QP) algorithm. See [1](#references) for details.

### Set up

SignatureEstimationPy requires Python 3. We recommend using Anaconda to install dependencies.

Once you have Anaconda installed, you can create an environment with all the dependencies installed, and activate it, with the following commands:

    conda env create -f environment.yml
    source activate signature-estimation-py-env

## Usage

The input to the `signature_estimation.py` script is a mutation count table and a mutation signatures table, and an (optional) list of signatures. To get the full list of parameters, run:

   python signature_estimation.py -h

## Example

We provide scripts in [`example/`](https://github.com/lrgr/signature-estimation-py/tree/master/example) to run SignatureEstimation on the [COSMIC mutation signatures](http://cancer.sanger.ac.uk/cosmic/signatures) and mutation data from the [Alexandrov, et al. (Nature 2013)](https://www.nature.com/articles/nature12477) paper.

Within the `example/` directory, you can run a full example on the 397 melanoma samples from the Alexandrov, et al. data using the provided `Snakefile` with the command:

    snakemake all

This will download and process the COSMIC signatures and melanoma mutation data, and compute their exposures using the SignatureEstimation quadratic program. The pipeline will produce a plot of exposures per signature, as automatically produced via continuous-integration below.

<img src='https://signature-estimation-py.lrgr.io/Alexandrov-et-al-Melanoma-exomes-cosmic-exposures.svg'>

## Acknowledgements

We gratefully acknowledge Xiaoqing Huang for her help developing and testing this code.

## References

Xiaoqing Huang, Damian Wojttowicz, and Teresa M Przytycka. Detecting presence of mutational signatures in cancer with confidence. *Bioinformatics(Oxford, England)*, 2017. doi: [10.1093/bioinformatics/btx604](https://doi.org/10.1093/bioinformatics/btx604).
