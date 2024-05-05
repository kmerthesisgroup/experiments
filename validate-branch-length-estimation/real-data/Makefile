DATASET := dataset
RESULT := result

$(RESULT): $(DATASET)
	docker run -v $(CURDIR):/workspace -w /workspace --rm kmerthssgrp/kmerfreqbdps python3 run.py

$(DATASET):
	mkdir -p $(DATASET)
	docker run --rm -v $(CURDIR)/$(DATASET):/dataset -v $(CURDIR)/gen-dataset.sh:/gen-dataset.sh kmerthssgrp/generate-kmer-count-matrix /gen-dataset.sh
