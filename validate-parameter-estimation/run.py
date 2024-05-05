#!/usr/env/python

import os
import yaml

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import phylogenetic_tree as pt
import bdps_parameter_estimator as bpt


np.random.seed(0)

RESULTS_DIR = 'result'
DATASET_DIR = 'dataset'
LAMBDA_BY_MU_IMAGE = f'{RESULTS_DIR}/lambda-by-mu-estimated-vs-original'
M_BY_MU_IMAGE = f'{RESULTS_DIR}/m-by-mu-estimated-vs-original'
CORREF_FILE =  f'{RESULTS_DIR}/corref.yaml'


def get_datasets(dataset_dir):
    data_files = []
    for root, dirs, files in os.walk(dataset_dir):
        for file in files:
            file_path = os.path.join(root, file)
            data_files.append(file_path)

    return data_files


def process_tree(filename,only_leaves):
    tree = pt.load_from_file(filename,only_leaves)
    ds = bpt.get_dataset(tree)
    estimator = bpt.BDPSParameterEstimator(ds)
    a,b,c = (estimator.estimate_parameter(bpt.initializer),
             ds.lamda/ds.mu,ds.m/ds.mu)

    el = np.abs((a[0]-b)/b)*100
    em = np.abs((a[1]-c)/c)*100
    e = np.sqrt(el**2+em**2)

    theta = [a[0],a[1]/a[0]]
    likelihood = -estimator.likelihood_func(theta)
    likelihood_real = -estimator.likelihood_func([ds.lamda/ds.mu,ds.m/ds.lamda])

    return {
                'filename':filename,
                'lamda':ds.lamda,
                'mu':ds.mu,
                'm':ds.m,
                'coeff_{lamda/mu}':a[0],
                'coeff_{m/mu}':a[1],
                'real_likelihood':likelihood_real,
                'estimated_likelihood':likelihood,
                'error':e
            }


def draw_lambda_by_mu_curve(df):
    lambda_orig = df['lamda']/df['mu']
    lambda_est = df['coeff_{lamda/mu}']

    plt.figure(figsize=(8, 8))

    plt.scatter(x=lambda_orig,y=lambda_est,marker='x')
    plt.plot([0,1],[0,1],color='r')

    plt.xlabel("original values of $\dfrac{\lambda}{\mu}$",fontsize=18)
    plt.ylabel("estimated values of $\dfrac{\lambda}{\mu}$",fontsize=18)

    plt.title("estimated vs original values of $\dfrac{\lambda}{\mu}$",fontsize=20)

    print('Saving lambda/mu image')
    plt.savefig(f'{LAMBDA_BY_MU_IMAGE}.png')
    plt.savefig(f'{LAMBDA_BY_MU_IMAGE}.eps', format='eps')
    print('Done')


def draw_m_by_mu_curve(df):
    #removing outlier (result still stands, it's just is really big and makes plot lose details)
    df = df[df['m']/df['mu'] < 20]
    m_orig = df['m']/df['mu']
    m_est = df['coeff_{m/mu}']


    plt.figure(figsize=(8, 8))

    plt.scatter(x=m_orig,y=m_est,marker='x')
    plt.plot([0,18],[0,18],color='r')

    plt.xlabel("original values of $\dfrac{m}{\mu}$",fontsize=18)
    plt.ylabel("estimated values of $\dfrac{m}{\mu}$",fontsize=18)


    plt.title("estimated vs original values of $\dfrac{m}{\mu}$",fontsize=20)


    print('Saving m/mu image')
    plt.savefig(f'{M_BY_MU_IMAGE}.png')
    plt.savefig(f'{M_BY_MU_IMAGE}.eps', format='eps')
    print('Done')


def generate_correlation_coefficients_file(df):
    results = {}

    lambda_orig = df['lamda']/df['mu']
    lambda_est = df['coeff_{lamda/mu}']
    results['corref-between-lambda/mu-estimated-and-original'] = float(np.corrcoef(lambda_est, lambda_orig)[0,1])

    m_orig = df['m']/df['mu']
    m_est = df['coeff_{m/mu}']
    results['corref-between-m/mu-estimated-and-original'] = float(np.corrcoef(m_est, m_orig)[0,1])

    with open(CORREF_FILE, 'w') as file:
        yaml.dump(results, file, default_flow_style=False)

    print('Correlation Coefficients:')
    print(yaml.dump(results, default_flow_style=False))


def main():
    tree_files = get_datasets(DATASET_DIR)

    df_parameters = pd.DataFrame(columns=[
        'filename',
        'lamda','mu','m','coeff_{lamda/mu}',
        'coeff_{m/mu}','real_likelihood',
        'estimated_likelihood','error'
    ])
    df_parameters.index.name = 'No'

    for file in tree_files:
        print(file, flush=True)

        entry = process_tree(file,only_leaves=False)
        df_parameters = df_parameters.append(entry, ignore_index=True)
        print(entry, flush=True)



    draw_lambda_by_mu_curve(df_parameters)
    draw_m_by_mu_curve(df_parameters)
    generate_correlation_coefficients_file(df_parameters)

if __name__ == '__main__':
    #generate results dir
    try:
        os.makedirs(RESULTS_DIR)
    except OSError as e:
        if not os.path.isdir(RESULTS_DIR):
            raise

    main()
