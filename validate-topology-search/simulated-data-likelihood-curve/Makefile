DATASET := dataset
INPUT_TOPOLOGY := input-topology
RESULT := result
IMAGE := simulated-data-likelihood-curve.eps

$(IMAGE): $(RESULT)
	docker run -v $(CURDIR):/workspace -w /workspace --rm kmerthssgrp/kmerfreqbdps python3 draw-curve.py

$(RESULT): $(DATASET)
	docker run -v $(CURDIR):/workspace -w /workspace --rm kmerthssgrp/kmerfreqbdps python3 compute-likelihood.py

$(DATASET): $(INPUT_TOPOLOGY)
	mkdir -p $(DATASET)
	docker run --rm -v $(CURDIR)/$(DATASET):/dataset -v $(CURDIR)/$(INPUT_TOPOLOGY):/input-topology -v $(CURDIR)/generate.py:/generate.py kmerthssgrp/kmerfreqbdps python3 /generate.py

$(INPUT_TOPOLOGY):
	docker run -v $(CURDIR):/workspace -w /workspace --rm ubuntu:jammy-20230804 bash -c "apt update && apt install -y curl make && cd generate-input-topology && make"
