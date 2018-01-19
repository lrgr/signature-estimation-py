configfile: 'config.yml'
from os.path import join
from urllib.parse import quote

################################################################################
# SETTINGS, FILES, AND DIRECTORIES
################################################################################
# Settings
SEQ_TYPE = config['seq_type'] + 's'
QUOTED_CANCER_TYPE = quote(config['cancer_type'])
CANCER_TYPE = config['cancer_type'].replace(' ', '')

# Directories
DATA_DIR = 'data'
SIGNATURES_DIR = join(DATA_DIR, 'signatures')
MUTATIONS_DIR = join(DATA_DIR, 'mutations')
RAW_DIR = join(DATA_DIR, 'raw')

# Files
RAW_COSMIC_SIGNATURES = join(RAW_DIR, 'cosmic_signatures_probabilities.txt')
COSMIC_SIGNATURES_JSON = join(SIGNATURES_DIR, 'cosmic-signatures.json')

RAW_MUTATION_CATALOGUE = join(RAW_DIR, '%s_%s_mutational_catalog_96_subs.txt' % (CANCER_TYPE, SEQ_TYPE))
PANDAS_MUTATION_CATALOGUE = join(MUTATIONS_DIR, '%s.tsv' % CANCER_TYPE)

COSMIC_ALEXANDROV_EXAMPLE = '%s-alexandrov-et-al-cosmic-exposusures.tsv' % config['cancer_type'].lower()

################################################################################
# RULES
################################################################################
# Data
rule download_cosmic_signatures:
    params:
        url='http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt'
    output:
        RAW_COSMIC_SIGNATURES
    shell:
        'wget -O {output} {params.url}'

rule process_cosmic_signatures:
    input:
        RAW_COSMIC_SIGNATURES
    output:
        COSMIC_SIGNATURES_JSON
    shell:
        'python data/process_cosmic_signatures.py -i {input} -o {output}'

rule download_mutations:
    params:
        url='ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/mutational_catalogs/%s/%s/%s_%s_mutational_catalog_96_subs.txt' % (SEQ_TYPE, QUOTED_CANCER_TYPE, QUOTED_CANCER_TYPE, SEQ_TYPE)
    output:
        RAW_MUTATION_CATALOGUE
    shell:
        'wget -O {output} {params.url}'

rule process_mutations:
    input:
        RAW_MUTATION_CATALOGUE
    output:
        PANDAS_MUTATION_CATALOGUE
    shell:
        'python data/process_alexandrov_et_al_mutations.py -i {input} -o {output}'


# Example run
rule example:
    input:
        signatures=COSMIC_SIGNATURES_JSON,
        mutations=PANDAS_MUTATION_CATALOGUE
    params:
        signatures='' if config['signatures'] is None else '-s ' + ' '.join([ '"%s"' % s for s in config['signatures']])
    output:
        COSMIC_ALEXANDROV_EXAMPLE
    shell:
        'python signature_estimation_qp.py -mf {input.mutations} '\
        '-sf {input.signatures} {params.signatures} -o {output}'
# All
rule all:
    input:
        COSMIC_ALEXANDROV_EXAMPLE
