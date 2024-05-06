DATASET := dataset
INPUT_TOPOLOGY := input-topology
RESULT := result
IMAGE := seven-primates-likelihood-curve.eps

$(IMAGE): $(RESULT)
	docker run -v $(CURDIR):/workspace -w /workspace --rm kmerthssgrp/kmerfreqbdps python3 draw-curve.py

$(RESULT): $(DATASET)
	docker run -v $(CURDIR):/workspace -w /workspace --rm kmerthssgrp/kmerfreqbdps python3 compute-likelihood.py

$(DATASET): $(INPUT_TOPOLOGY)
	mkdir -p $(DATASET)
	docker run --rm -v $(CURDIR)/$(DATASET):/dataset -v $(CURDIR)/$(INPUT_TOPOLOGY):/input-topology -v $(CURDIR)/generate.py:/generate.py kmerthssgrp/kmerfreqbdps python3 /generate.py

$(INPUT_TOPOLOGY):
	mkdir -p $(INPUT_TOPOLOGY)
	docker run --rm -v $(CURDIR)/$(INPUT_TOPOLOGY):/input-topology -v $(CURDIR)/generate-input-topology.sh:/generate-input-topology.sh kmerthssgrp/generate-kmer-count-matrix /generate-input-topology.sh
