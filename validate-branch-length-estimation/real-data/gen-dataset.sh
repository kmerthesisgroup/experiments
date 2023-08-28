#!/bin/bash
set -xe

cd seven-primates
cp ../run.sh ../generate-matrix.sh ../clean.sh ./

LOWER=7
HIGHER=13
NUM_SAMPLE=2000

DEST_DIR=/dataset

for (( i=$LOWER; i<=$HIGHER; i++ ))
do
	./generate-matrix.sh seven-primates $i $NUM_SAMPLE
	./run.sh seven-primates $i $NUM_SAMPLE
	cp original_topology-seven-primates-$i-$NUM_SAMPLE/final_tree.txt $DEST_DIR/seven-primates-$i-$NUM_SAMPLE-tree.txt
	cp original_topology-seven-primates-$i-$NUM_SAMPLE/final_species_file.txt $DEST_DIR/seven-primates-$i-$NUM_SAMPLE-species.txt
	cp seven-primates-$i-$NUM_SAMPLE/entropy.txt $DEST_DIR/seven-primates-$i-$NUM_SAMPLE-entropy.txt
done
