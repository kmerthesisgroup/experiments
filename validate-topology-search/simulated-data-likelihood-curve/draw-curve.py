import os
import statistics
import yaml

import matplotlib.pyplot as plt

RESULT_DIR = 'result'
IMAGE_NAME = 'simulated-dataset-likelihood-curve'

def load_data():
    files = os.listdir(RESULT_DIR)
    data = dict()
    for file in files:
        if file.endswith('-tree.yaml'):
            splt = file[:-len('-tree.yaml')].split('-')
            rf_distance = int(splt[-2])
            if rf_distance not in data:
                data[rf_distance] = list()

            full_filepath = os.path.join(RESULT_DIR, file)
            with open(full_filepath, 'r') as fp:
                yaml_data = yaml.safe_load(fp)
                data[rf_distance].append(yaml_data['estimated_tree_likelihood'])

    return data

def get_mean_and_std(data):
    keys = sorted(data)
    mu = [ statistics.mean(data[key]) for key in keys]
    std = [ statistics.stdev(data[key]) for key in keys]

    return keys, mu, std

def draw_image(keys, mu, std):
    plt.figure(figsize=(8, 6))
    plt.plot(keys, mu, 'b-', label='Line')

    plt.errorbar(keys, mu, yerr=std, fmt='none', ecolor='gray', alpha=0.4, capsize=3, label='Error Bars')
    plt.legend(fontsize=18)

    plt.xlabel('RF Distance', fontsize=18)
    plt.ylabel('Likelihood' , fontsize=18)
    plt.title('Likelihood vs RF Distance' , fontsize=20)

    plt.savefig(f'{IMAGE_NAME}.png')
    plt.savefig(f'{IMAGE_NAME}.eps', format='eps')

def main():
    data = load_data()
    key, mu, std = get_mean_and_std(data)
    draw_image(key, mu, std)

if __name__ == '__main__':
    main()
