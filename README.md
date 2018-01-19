# Signature Estimation in Python

Adapting the [Signature Estimation R](https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/index.cgi#signatureestimation) package for Python.

### Set up

SignatureEstimationPy requires Python 3. We recommend using Anaconda to install dependencies.

Once you have Anaconda installed, you can create an environment with all the dependencies installed, and activate it, with the following commands:

    conda env create -f environment.yml
    source activate signature-estimation-py-env

### Usage

The input to the `signature_estimation_qp.py` script is a mutation count table and a JSON file of mutation signatures, and an (optional) list of signatures.

#### Example

We provide example scripts for downloading and processing the [COSMIC mutation signatures](http://cancer.sanger.ac.uk/cosmic/signatures) and mutation data from the [Alexandrov, et al. (Nature 2013)](https://www.nature.com/articles/nature12477) paper in [`data/`](https://github.com/lrgr/signature-estimation-py/tree/master/data).

You can run a full example on the 397 melanoma samples from the Alexandrov, et al. data using the provided `Snakefile` with the command:

    snakemake all

This will download and process the COSMIC signatures and melanoma mutation data, and compute their exposures using the SignatureEstimation quadratic program.

#### Configuration

You can configure the `Snakefile` by changing the values of parameters in `config.yml`, e.g. to run on whole-genomes or other cancer types. When the `signatures` parameter is removed entirely, `signature_estimation_qp.py` will compute exposures for all provided signatures.
