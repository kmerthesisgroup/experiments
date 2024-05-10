import ete3
import os
import argparse
import yaml

def get_distance(newick1, newick2):
    t1 = ete3.Tree(newick1)
    t2 = ete3.Tree(newick2)

    dist = t1.robinson_foulds(t2, unrooted_trees=True)[0]

    return dist

def get_trees(root_dir):
    files = os.listdir(root_dir)

    trees = dict()
    for file in files:
        full_path = os.path.join(root_dir, file)
        with open(full_path, 'r') as fp:
            tree = fp.read()
            trees[file] = tree

    return trees

def cluster_trees(trees):
    clusters = []
    for filename, tree in trees.items():
        done = False
        for cluster in clusters:
            t1 = trees[cluster[0]]
            dist = get_distance(t1, tree)
            if dist == 0:
                cluster.append(filename)
                done = True
                break
        if not done:
            clusters.append([filename])

    return clusters

def main(root_dir):
    trees = get_trees(root_dir)
    clusters = cluster_trees(trees)
    yaml_data = { f'cluster-{i}' : clusters[i] for i in range(len(clusters))}

    outfile = os.path.basename(root_dir)
    if outfile == '':
        outfile = os.path.basename(root_dir[:-1])
    if outfile == '':
        outfile = 'out.yaml'
    else:
        outfile = outfile + '-clusters.yaml'

    with open(outfile, 'w') as file:
        yaml.dump(yaml_data, file, default_flow_style=False)

    tree_names = trees.keys()
    for tree1 in tree_names:
        for tree2 in tree_names:
            print(tree1, tree2, get_distance(trees[tree1], trees[tree2]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('root_dir')

    args = parser.parse_args()
    main(args.root_dir)
