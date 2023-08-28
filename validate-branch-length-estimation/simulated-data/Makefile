DATASET := dataset

.PHONY: clean

result: $(DATASET) 
	docker run -v $(CURDIR):/workspace -w /workspace --rm kmerthssgrp/kmerfreqbdps python3 run.py
		
$(DATASET):  
	docker run -v $(CURDIR):/workspace -w /workspace --rm ubuntu:jammy-20230804 bash -c "apt update && apt install -y curl make && cd generate-dataset && make" 

clean:
	cd generate-dataset && make clean
	rm -rf result
	rm -rf log 
