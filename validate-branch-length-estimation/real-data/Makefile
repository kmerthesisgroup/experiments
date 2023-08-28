DATASET := dataset
RESULT := result

$(RESULT): $(DATASET) 
	docker run -v $(CURDIR):/workspace -w /workspace --rm kmerthssgrp/kmerfreqbdps python3 run.py

$(DATASET): 
	mkdir -p $(DATASET)
	docker run --rm -v ./$(DATASET):/dataset -v ./gen-dataset.sh:/gen-dataset.sh kmerthssgrp/generate-kmer-count-matrix /gen-dataset.sh

