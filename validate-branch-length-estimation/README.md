# Description
- Generates datasets
  - in the simulated case, generates random trees
  - in case of real case, the seven primates dataset is used to generate a tree file that complies with our format. 
- Randomize the branch lengths of the tree and then use our model to optimize branch lengths
- Draw plots of estimated and actual branch lengths (in the same scale)
- Generate yaml file with likelihood and bldk
- Entropy file will be present in the dataset.

# How to Run
```
make dataset
make result
```

# Artifacts
- `tree.yaml`
- `tree-estimated.eps`
- Optimized tree (at `result/{filename}-{final-iteration}.tree`)
  
