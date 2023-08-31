#!/bin/bash
set -xe

DATASET=dataset
GENERATE_RANDOM_TREE=generate-random-tree-linux-amd64
RANDOM_TREE_CONF=random-tree.conf
RANDOM_TREE=$DATASET/random_topology-tree.txt

ORIGINAL_TREE=$(ls $DATASET/original_topology*)
LOWER=0.1
UPPER=10

SITES=$(head -n 2 $ORIGINAL_TREE | tail -n 1)
SPECIES=$(head -n 1 $ORIGINAL_TREE)

cat << EOF > $RANDOM_TREE_CONF
seed: 0

branch-length-distribution:
  uniform:
    lower: 0.1
    upper: 1

lambda: -1
mu: -1
m: -1

number-of-species: $SPECIES
number-of-sites: $SITES
EOF

if [ ! -f $GENERATE_RANDOM_TREE ]
then
	curl -LJO https://github.com/hsiam261/rand-phylo-tree/releases/download/v0.0.1/generate-random-tree-linux-amd64
fi
chmod +x $GENERATE_RANDOM_TREE
./$GENERATE_RANDOM_TREE $RANDOM_TREE_CONF $RANDOM_TREE

random_string=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 10 | head -n 1)
head -n $(( 2*SPECIES + 1 )) $RANDOM_TREE > /tmp/${random_string}.txt
mv /tmp/${random_string}.txt $RANDOM_TREE

tail -n +$(( 2*SPECIES + 2 )) $ORIGINAL_TREE >> $RANDOM_TREE
