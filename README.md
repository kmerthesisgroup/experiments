a collection of experiments conducted during writing of the paper "A Birth-Death-Migration Model for Alignment-Free Phylogeny Estimation using k-mer Frequencies" 

## Experiments
- validate parameter estimation
- validate branch length estimation
  - simulated data
  - real data   
- validate topology search
  - simulated data
  - real data 
- afproject experiments
  - ecoli
  - plant
  - fish
 
## Setup
The might require the following docker images locally present in the machine:
- `kmerthssgrp/kmerfreqbdps`
- `kmerthssgrp/generate-kmer-count-matrix`
Since they aren't published in any public registry, they need to be built first.
The experiments will also require bash, make and docker to be installed on the machine.

So best way to set up environment for all experiments is to do the following:
- create a folder called `workspace`
- cd into the directory and clone the following three repositories:
  - https://github.com/kmerthesisgroup/generate-kmer-count-matrix
  - https://github.com/kmerthesisgroup/experiments
  - https://github.com/kmerthesisgroup/kwordfreqbdps
- build the docker images

We can do this in ubuntu using the following commands
```
sudo apt update
sudo apt install make git

# Install Docker
# Add Docker's official GPG key:
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update

sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

# create workspace
mkdir workspace
cd workspace

git clone git@github.com:kmerthesisgroup/experiments.git
git clone git@github.com:kmerthesisgroup/kwordfreqbdps.git
git clone git@github.com:kmerthesisgroup/generate-kmer-count-matrix.git

# prepare docker images
cd generate-kmer-count-matrix
docker build -t kmerthssgrp/generate-kmer-count-matrix .
cd ../

cd kmerfreqbdps
docker build -t kmerthssgrp/kmerfreqbdps .
cd ../
```
