#!/usr/bin/env python3
"""Usage:
  fit-inferelator.py <prefix> <method> <x> <reps>

Arguments:
  <prefix>  Input data is <prefix>-exprs.tsv, <prefix>-prior.tsv, ...
  <method>  bbsr or enet
  <x>       prior strength (bbsr) or alpha (enet)
  <reps>    Replications (seeds)
"""

from docopt import docopt
arg = docopt(__doc__)
from inferelator import inferelator_workflow, inferelator_verbose_level, MPControl, CrossValidationManager

# Set verbosity level to "Talky"
inferelator_verbose_level(1)

# Set the location of the input data and the desired location of the output files
DATA_DIR = 'network/inferelator'
OUTPUT_DIR = 'network/inferelator/'
EXPRESSION_FILE_NAME = arg['<prefix>'] + '-exprs.tsv'
PRIORS_FILE_NAME = arg['<prefix>'] + '-prior.tsv'
TF_LIST_FILE_NAME = arg['<prefix>'] + '-tf.txt'

CV_SEED = list(range(1, int(arg['<reps>'])+1))

import os

n_cores_local = int(os.environ.get("NSLOTS", os.cpu_count()))

if __name__ == '__main__':
    MPControl.set_multiprocess_engine("multiprocessing", processes=n_cores_local)
    MPControl.connect()

def set_up_workflow(wkf):
    wkf.set_file_paths(
        input_dir=DATA_DIR,
        output_dir=OUTPUT_DIR,
        tf_names_file=TF_LIST_FILE_NAME,
        priors_file=PRIORS_FILE_NAME
    )
    wkf.set_expression_file(tsv=EXPRESSION_FILE_NAME)
    wkf.set_file_properties(expression_matrix_columns_are_genes=False)
    wkf.set_network_data_flags(use_no_gold_standard=True)
    wkf.set_run_parameters(num_bootstraps=5)
    wkf.set_crossvalidation_parameters(
        split_gold_standard_for_crossvalidation=True,
        cv_split_ratio=0.2
    )
    return wkf

# Create a worker
if arg['<method>'] == 'bbsr':
    worker = inferelator_workflow(regression="bbsr", workflow="single-cell")
    worker = set_up_workflow(worker)
    worker.append_to_path("output_dir", arg['<prefix>'])
    worker.set_regression_parameters(
        prior_weight=1 + float(arg['<x>'])
    )
elif arg['<method>'] == 'enet':
    worker = inferelator_workflow(regression="elasticnet", workflow="single-cell")
    worker = set_up_workflow(worker)
    worker.append_to_path("output_dir", arg['<prefix>'])
    worker.set_regression_parameters(
        l1_ratio=float(arg['<x>']),
        max_iter=2000,
        eps=0.1
    )

# Create a crossvalidation wrapper
cv_wrap = CrossValidationManager(worker)

# Assign variables for grid search
cv_wrap.add_gridsearch_parameter('random_seed', CV_SEED)

# Run
cv_wrap.run()
