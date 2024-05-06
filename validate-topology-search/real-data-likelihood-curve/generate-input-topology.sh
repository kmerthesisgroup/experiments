#!/bin/bash
set -xe

cd seven-primates
cp ../run.sh ../generate-matrix.sh ../clean.sh ./

NUM_SAMPLE=2000
KMER_LENGTH=8 #yeilded the highest entropy

DEST_DIR=/input-topology

./generate-matrix.sh seven-primates $KMER_LENGTH $NUM_SAMPLE
./run.sh seven-primates $KMER_LENGTH $NUM_SAMPLE
cp original_topology-seven-primates-$KMER_LENGTH-$NUM_SAMPLE/final_tree.txt $DEST_DIR/seven-primates-$KMER_LENGTH-$NUM_SAMPLE-tree.txt
cp original_topology-seven-primates-$KMER_LENGTH-$NUM_SAMPLE/final_species_file.txt $DEST_DIR/seven-primates-$KMER_LENGTH-$NUM_SAMPLE-species.txt
cp seven-primates-$KMER_LENGTH-$NUM_SAMPLE/entropy.txt $DEST_DIR/seven-primates-$KMER_LENGTH-$NUM_SAMPLE-entropy.txt
