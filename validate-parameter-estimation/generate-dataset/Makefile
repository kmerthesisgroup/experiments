GENERATE_RANDOM_TREE := ./generate-random-tree-linux-amd64 

DATASET := ../dataset
CONF_DIR := conf

# List all dataset configurations from source directory
DATASET_CONFS := $(wildcard $(CONF_DIR)/*.yaml)

# All dataset directories
DATASET_DIRS := $(patsubst $(CONF_DIR)/%.yaml, $(DATASET)/%/, $(DATASET_CONFS))

.PHONY: generate_data clean

generate-random-tree-linux-amd64:
	curl -LJO https://github.com/hsiam261/rand-phylo-tree/releases/download/v0.0.1/generate-random-tree-linux-amd64
	chmod +x $(GENERATE_RANDOM_TREE)

$(DATASET_DIRS): $(DATASET)/%/: $(CONF_DIR)/%.yaml generate-random-tree-linux-amd64
	mkdir -p $@
	$(GENERATE_RANDOM_TREE) $< $@tree
 
generate_data: $(DATASET_DIRS)

clean:
	rm -rf $(DATASET)
	rm -f generate-random-tree-linux-amd64
