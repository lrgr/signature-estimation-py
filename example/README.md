# Example

Scripts to run SignatureEstimation on COSMIC mutation signature and data from the Alexandrov, et al. (Nature, 2013) paper. To run the full example, execute:

    snakemake all

## Configuration

By default, we just run on the Alexandrov, et al. melanoma exomes. You can [configure the `Snakefile`](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html) to run on whole-genomes or -exomes (parameter `seq_type`), other cancer types (`cancer_type`), or to run on specific signatures (`signatures`). By default, all COSMIC signatures are used.

