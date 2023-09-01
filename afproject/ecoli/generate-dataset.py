#!/usr/env/python

import subprocess
import argparse
import os

FASTA_FILE_DIR='data'

SPECIES='fish'
LOWER_K=5
HIGHER_K=17

LOWER_NUM_SAMPLE=5000
HIGHER_NUM_SAMPLE=50000

DEST_DIR='dest'


def find_highest_num_sample():
    longest_size = 0
    for fasta_file in os.listdir(FASTA_FILE_DIR):
        file_path = os.path.join(FASTA_FILE_DIR, fasta_file)
        if os.path.isfile(file_path) and file_path.endswith('.fasta'):
            file_size = os.path.getsize(file_path)
            if(longest_size < file_size):
                longest_size = file_size
    print(longest_size)
    return longest_size


def optimize_hyper_parameters():
    best_entropy = -1
    best_k = -1
    best_num_sample = -1
    entropies = []
    for k in range(LOWER_K, HIGHER_K+1):
        num_sample = int(min(max((4**(k-1)/k)//1000*1000, LOWER_NUM_SAMPLE), HIGHER_NUM_SAMPLE))
        print(subprocess.check_output(f'./generate-matrix.sh {SPECIES} {k} {num_sample}', shell=True).decode('utf-8'))	
        entropy=float(subprocess.check_output(f'cat {SPECIES}-{k}-{num_sample}/entropy.txt', shell=True).decode('utf-8'))/num_sample
        entropies.append((k, entropy, num_sample))
        if entropy >= best_entropy:
            best_entropy = entropy
            best_k = k
            best_num_sample = num_sample
    
    with open('entropy.txt', 'w') as fp:
        for e,k,n in entropies:
            print(e,k,n, file=fp)

    return (best_entropy, best_k, best_num_sample)

def main():
    subprocess.check_output(f'./download-{SPECIES}.sh', shell=True)
    subprocess.check_output('cp ../run.sh ../generate-matrix.sh ../clean.sh ./', shell=True)

    global HIGHER_NUM_SAMPLE
    HIGHER_NUM_SAMPLE=min(HIGHER_NUM_SAMPLE, find_highest_num_sample())
    e, k, n = optimize_hyper_parameters()

    print("Highest Entropy/num_sample : ", e)
    print("k : ", k)
    print("num_samples : ", n)
    print(subprocess.check_output(f'./run.sh {SPECIES} {k} {n}', shell=True).decode('utf-8'))

    subprocess.check_output(f'mkdir -p {DEST_DIR}', shell=True)
    
    for file in os.listdir('newick'):
        base_name = file.split('.')[0]
        subprocess.check_output(f'cp {base_name}-{SPECIES}-{k}-{n}/final_tree.txt {DEST_DIR}/{base_name}-{k}-{n}-tree.txt', shell=True)
    
    subprocess.check_output(f'cp {SPECIES}-{k}-{n}/species-file.txt {DEST_DIR}/', shell=True)
    subprocess.check_output(f'cp entropy.txt {DEST_DIR}/', shell=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("--species", action="store", required=True,type=str, help="name of species")
    parser.add_argument("--klower", action="store", required=True, type=int, help="lowest value of k")
    parser.add_argument("--khigher", action="store", required=True, type=int, help="highest value of k")
    parser.add_argument("--sample-lower", action="store", required=True, type=int, help="lowest value of num-sample")
    parser.add_argument("--sample-higher", action="store", required=True, type=int, help="highest valueof num-sample")

    args = parser.parse_args()
    SPECIES = args.species
    LOWER_K = args.klower
    HIGHER_K = args.khigher
    LOWER_NUM_SAMPLE = args.sample_lower
    HIGHER_NUM_SAMPLE = args.sample_higher

    main()
