# Signature Estimation in Python

Adapting the [Signature Estimation R](https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/index.cgi#signatureestimation) package for Python.

### Set up

SignatureEstimationPy requires Python 3. We recommend using Anaconda to install dependencies.

Once you have Anaconda installed, you can create an environment with all the dependencies installed, and activate it, with the following commands:

    conda env create -f environment.yml
    source activate signature-estimation-py-env

### Usage

The input to the `signature_estimation_qp.py` script is a mutation count table and a JSON file of mutation signatures. TO-DO: More documentation on the format of these files.
