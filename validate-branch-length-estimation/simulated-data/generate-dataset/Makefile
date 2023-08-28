GENERATE_RANDOM_TREE := ./generate-random-tree-linux-amd64

DATASET := ../dataset
CONF := 10-species-2000-sites-small-branches-normalized-params.yaml
RESULT_DIR := result

.PHONY: default clean

default: $(DATASET)

$(GENERATE_RANDOM_TREE):
	curl -LJO https://github.com/hsiam261/rand-phylo-tree/releases/download/v0.0.1/generate-random-tree-linux-amd64
	chmod +x $(GENERATE_RANDOM_TREE)

$(DATASET): $(GENERATE_RANDOM_TREE) $(CONF)
	mkdir -p $@
	$(GENERATE_RANDOM_TREE) $(CONF) $(DATASET)/tree
 
clean:
	rm -rf $(DATASET)
	rm -f generate-random-tree-linux-amd64
