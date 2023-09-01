#!/usr/env/bin python

import os
import yaml
import io
import random
import argparse
import numpy as np
import matplotlib.pyplot as plt

import phylogenetic_tree as pt
import felsenstein_pruner as fpruner
import branch_length_estimator as ble
import bdps_parameter_estimator as bpe

from functools import lru_cache
from Bio import Phylo

np.random.seed(0)
random.seed(0)

LOG_DIR = 'log'
RESULT_DIR = 'result'
DATASET_DIR = 'dataset'

LOWER=0.001
UPPER=20
NUM_PASSES=30
EPRECISSION=np.float32(0.001)
LPRECISSION=np.float32(0.01)
MIN_MAX_KMER_COUNT=20
MAX_MAX_KMER_COUNT=500
CUT_OFF_PI=10**-8

@lru_cache
def get_max_kmer_count(tree):
    Q = fpruner.gen_censored_linear_bdps_qmat(tree.lamda, tree.m, tree.mu, MAX_MAX_KMER_COUNT)
    lnpi = bpe.gen_ln_pi_array(tree.lamda, tree.m/tree.lamda, MAX_MAX_KMER_COUNT)
    pi = np.e**lnpi
    pi = pi/np.sum(pi)
    return np.argwhere(pi>CUT_OFF_PI)[-1, 0]

@lru_cache
def estimate_parameters(tree):
    ds = bpe.get_dataset(tree)
    estimator = bpe.BDPSParameterEstimator(ds)
    a = estimator.estimate_parameter(bpe.initializer,number_of_trials=10)
    theta = [a[0],a[1]/a[0]]
    tree.lamda = a[0]
    tree.mu = 1
    tree.m = a[1]
    print('estimated_paramters: ', theta)

def get_likelihood(tree):
    Q = fpruner.gen_censored_linear_bdps_qmat(tree.lamda, tree.m, tree.mu, get_max_kmer_count(tree))
    lnpi = bpe.gen_ln_pi_array(tree.lamda, tree.m/tree.lamda, get_max_kmer_count(tree))
    pi = np.e**lnpi
    pi = pi/np.sum(pi)
    print(pi, flush=True)
    pruner = fpruner.FelsensteinPruner(tree, qmat=Q, pi=pi)
    return pruner.compute_log_likelihood(0, recompute_table=True)

def draw_trees(tree, filename):
    fig , ax = plt.subplots(figsize=(8, 6))

    pt.draw_tree(tree, ax)
    ax.set_title(f'Estimated Tree for {filename}', fontsize=20)
    ax.set_xlabel('branch length', fontsize=16)
    ax.set_ylabel('taxa', fontsize=16)

    plt.savefig(f'{RESULT_DIR}/{filename}.png')

def process_tree(tree_original_file, tree_file, species_file):
    tree_original = pt.load_from_file(tree_original_file, only_leaves=True)
    tree_original.load_species_file(species_file)

    tree = pt.load_from_file(tree_file, only_leaves=True)
    tree.load_species_file(species_file)
    for i in range(len(tree.edges)):
        tree.update_branch_length(i,np.random.rand()*5)

    result = dict()
    
    estimate_parameters(tree)

    file_name = os.path.basename(tree_file).split('.')[0]
    max_kmer_count = get_max_kmer_count(tree)
    print('max kmer count: ', max_kmer_count) 
    ble.estimate_branch_length(tree, LOWER, UPPER, max_kmer_count, num_passes=NUM_PASSES, seed_tree=True, 
                               eprecisson=EPRECISSION, lprecisson=LPRECISSION,logfile=f'{LOG_DIR}/{file_name}')
    pt.dump_to_file(tree, f'{RESULT_DIR}/{file_name}')
    print('Done dumping file')

    result['estimated_tree_likelihood'] = float(get_likelihood(tree))
    result['rf_distance'] = pt.rf_distance(tree_original, tree) 
    
    print(result)
    with open(f'{RESULT_DIR}/{file_name}.yaml', 'w') as fp:
        yaml.dump(result, fp, default_flow_style=False)
    

    tree.load_species_file(species_file)
    draw_trees(tree, file_name)
    
    return result

if __name__ == '__main__':
    #generate results dir
    try:
        os.makedirs(RESULT_DIR)
    except OSError as e:
        if not os.path.isdir(RESULT_DIR):
            raise
    # generate log dir
    try:
        os.makedirs(LOG_DIR)
    except OSError as e:
        if not os.path.isdir(LOG_DIR):
            raise

    parser = argparse.ArgumentParser()
    parser.add_argument("--species-file", action="store", required=True,type=str, help="path to species-file")
    parser.add_argument("--original-tree-file", action="store", required=True,type=str, help="path to original-tree-file")
    parser.add_argument("--tree-file", action="store", required=True,type=str, help="path to tree-file")
    args = parser.parse_args()
    
    process_tree(args.original_tree_file, args.tree_file, args.species_file) 
