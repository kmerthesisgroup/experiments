SPECIES := fish
KLOWER := 5
KHIGHER := 17
SAMPLE_LOWER := 5000
SAMPLE_HIGHER := 50000

DATASET := dataset
RESULT := result

TREES := $(wildcard $(DATASET)/*-tree.txt)
SPECIES_FILE := $(DATASET)/species-file.txt

OPTIMIZED_TREES := $(patsubst $(DATASET)/%-tree.txt, $(RESULT)/%-tree.txt, $(TREES))
ORIGINAL_TOPOLOGY := $(wildcard $(DATASET)/original_topology-*-tree.txt)

.phony: trees cleanall

$(DATASET): 
	mkdir -p $(DATASET)
	mkdir -p fasta_dir
	docker run -v ./$(DATASET):/workspace/$(SPECIES)/dest -v ./generate-dataset.py:/generate-dataset.py -v ./fasta_dir:/workspace/$(SPECIES)/data  -w /workspace/$(SPECIES)/ --rm kmerthssgrp/generate-kmer-count-matrix python3 /generate-dataset.py --species $(SPECIES) --klower $(KLOWER) --khigher $(KHIGHER) --sample-lower $(SAMPLE_LOWER) --sample-higher $(SAMPLE_HIGHER)
	docker run -v $(CURDIR):/workspace -w /workspace --rm ubuntu:jammy-20230804 bash -c "apt update && apt install -y curl && ./generate-random-tree.sh"

trees: $(OPTIMIZED_TREES)

$(OPTIMIZED_TREES): $(RESULT)/%-tree.txt: $(DATASET)/%-tree.txt
	docker run -v $(CURDIR):/workspace -w /workspace --rm kmerthssgrp/kmerfreqbdps python3 generate-tree.py --species-file $(SPECIES_FILE) --original-tree-file $(ORIGINAL_TOPOLOGY) --tree-file $<


cleanall:
	rm -rf $(DATASET)
	rm -rf $(RESULT)
