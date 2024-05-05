#!/usr/env/bin python

import os
import yaml
import io
import random
import numpy as np
import matplotlib.pyplot as plt

import phylogenetic_tree as pt
import felsenstein_pruner as fpruner
import branch_length_estimator as ble
import bdps_parameter_estimator as bpe

from Bio import Phylo

np.random.seed(0)
random.seed(0)

LOG_DIR = 'log'
RESULT_DIR = 'result'
DATASET_DIR = 'dataset'


LOWER=0.001
UPPER=10
NUM_PASSES=30
EPRECISSION=np.float32(0.001)
LPRECISSION=np.float32(0.01)
MIN_MAX_KMER_COUNT=20
MAX_MAX_KMER_COUNT=40

def get_max_kmer_count(tree):
    return max(min(np.max(tree.kmer_count)*2, MAX_MAX_KMER_COUNT), MIN_MAX_KMER_COUNT)

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
    pruner = fpruner.FelsensteinPruner(tree, qmat=Q, pi=pi)
    return pruner.compute_log_likelihood(0, recompute_table=True)

def draw_trees(tree_original, tree, tree_filename):
    fig , axs = plt.subplots(2, 1, figsize=(8, 12), sharex=True,
                          sharey=True)
    pt.draw_tree(tree, axs[1])
    axs[1].set_title('Estimated Tree for Seven Primates Dataset', fontsize=20)
    axs[1].set_xlabel('branch length', fontsize=16)
    axs[1].set_ylabel('taxa', fontsize=16)

    pt.draw_tree(tree_original, axs[0])
    axs[0].set_title('Original Tree for Seven Primates Dataset',fontsize=20)
    axs[0].set_xlabel('branch length', fontsize=16)
    axs[0].set_ylabel('taxa', fontsize=16)

    file_name = os.path.basename(tree_filename).split('.')[0]

    plt.savefig(f'{RESULT_DIR}/{file_name}.eps', format='eps')

def process_tree(tree_file, species_file):
    tree_original = pt.load_from_file(tree_file, only_leaves=True)
    tree_original.load_species_file(species_file)

    tree = pt.load_from_file(tree_file, only_leaves=True)
    tree.load_species_file(species_file)
    for i in range(len(tree.edges)):
        tree.update_branch_length(i,np.random.rand()*5)

    result = dict()

    estimate_parameters(tree)
    result['random_tree_likelihood'] = float(get_likelihood(tree))
    result['random_tree_bldk'] = float(pt.bldk(tree_original, tree, pt.find_optimal_K(tree_original,tree)))

    file_name = os.path.basename(tree_file).split('.')[0]
    max_kmer_count = get_max_kmer_count(tree)
    print('max kmer count: ', max_kmer_count)
    ble.estimate_branch_length(tree, LOWER, UPPER, max_kmer_count, num_passes=NUM_PASSES, seed_tree=True,
                               eprecisson=EPRECISSION, lprecisson=LPRECISSION,logfile=f'{LOG_DIR}/{file_name}')

    result['estimated_tree_likelihood'] = float(get_likelihood(tree))
    result['estimated_tree_bldk'] = float(pt.bldk(tree_original, tree, pt.find_optimal_K(tree_original,tree)))
    draw_trees(tree_original, tree, tree_file)

    print(result)
    with open(f'{RESULT_DIR}/{file_name}.yaml', 'w') as fp:
        yaml.dump(result, fp, default_flow_style=False)

    return result

def get_datasets(dataset_dir):
    files = [ file for file in os.listdir(dataset_dir) if os.path.isfile(os.path.join(dataset_dir, file))]
    tree_files = [ file for file in files if file.endswith('-tree.txt') ]
    species_files = [ file.removesuffix('-tree.txt') + '-species.txt' for file in tree_files ]

    return zip(tree_files, species_files)


def main():
    datasets = get_datasets(DATASET_DIR)
    for tree_file, species_file in datasets:
        process_tree(f'{DATASET_DIR}/{tree_file}', f'{DATASET_DIR}/{species_file}')

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

    main()
