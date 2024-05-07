from copy import deepcopy
import random
import collections
import numpy as np
import phylogenetic_tree as pt
import os
import re
import numpy as np

random.seed(0)
np.random.seed(0)

INPUT_TOPOLOGY_DIRECTORY='/input-topology'
OUTPUT_DIRECTORY='/dataset'

ITERATIOINS = 100
MAX_SPRS = 10
MAX_RF_DISTANCE = 10
MAX_CNT = 15

def remove_edge_from_vertex_adjacency(tree,u,id_to_remove):
  tree = deepcopy(tree)
  for i in range(len(tree.adjList[u])):
    e = tree.adjList[u,i]
    if e is None or e.id==-1:
      break
    if e.id == id_to_remove:
      tree.degree[u]-=1
      tree.adjList[u,i] = tree.adjList[u,tree.degree[u]]
      tree.adjList[u,tree.degree[u]] = None
      break
  return tree

def delete_edge(tree,edge_id):
  e = tree.edges[edge_id]
  tree_copy = deepcopy(tree)
  tree_copy = remove_edge_from_vertex_adjacency(tree_copy,e.src,edge_id)
  tree_copy = remove_edge_from_vertex_adjacency(tree_copy,e.dest,edge_id)
  return tree_copy

def add_edge_into_index(tree,u,v,w,index):
  tree_copy = deepcopy(tree)

  tree_copy.add_edge(u,v,w)
  new_id = tree_copy.edges[-1].id

  tree_copy.edges[-1].id = index
  tree_copy.edges[index] = tree_copy.edges[-1]
  tree_copy.edges = tree_copy.edges[:-1]

  for i in range(len(tree_copy.adjList[u])):
    e = tree_copy.adjList[u,i]
    if e is None or e.id==-1:
      break
    if e.id == new_id:
      e.id = index
      break

  for i in range(len(tree_copy.adjList[v])):
    e = tree_copy.adjList[v,i]
    if e is None or e.id==-1:
      break
    if e.id == new_id:
      e.id = index
      break
  print(tree_copy.adjList)


  return tree_copy

# Never prune the root
# Will break invariants
def prune_vertex(tree,v,id_to_insert_into):
  edges_deleted = []
  tree_copy = deepcopy(tree)
  for e in tree_copy.adjList[v]:
    if e is None or e.id == -1:
      break
    tree_copy = delete_edge(tree_copy,e.id)
    edges_deleted.append(e.id)

  if len(edges_deleted) == 2:
    e = tree.edges[edges_deleted[0]]
    u = e.src
    if(u==v):
      u = e.dest

    e = tree.edges[edges_deleted[1]]
    p = e.src
    if(p==v):
      p = e.dest

    tree_copy = add_edge_into_index(tree_copy,u,p,1,id_to_insert_into)
    # tree_copy.add_edge(u,p,0)
    # tree_copy.edges[-1].id = id_to_insert_into
    # tree_copy.edges[id_to_insert_into] = tree.edges[-1]
    # tree_copy.edges = tree_copy.edges[:-1]

  return tree_copy,edges_deleted

def prune_subtree(tree,edge_id,endpoint_to_cut):
  tree_copy = deepcopy(tree)

  edge = tree_copy.edges[edge_id]


  other_endpoint = edge.src
  if other_endpoint == endpoint_to_cut:
    other_endpoint = edge.dest

  id_to_remove = edge.id
  tree_copy = delete_edge(tree_copy,id_to_remove)
  tree_copy, deleted_edges = prune_vertex(tree_copy,other_endpoint,id_to_remove)
  return tree_copy,other_endpoint,deleted_edges

def graft_subtree(tree,graft_edge_id,grafted_subtree_root, deleted_vertex, deleted_edge_ids):
  edge_to_graft_on = tree.edges[graft_edge_id]

  tree_copy = delete_edge(tree,graft_edge_id)
  tree_copy = add_edge_into_index(tree_copy,edge_to_graft_on.src,deleted_vertex,1,deleted_edge_ids[0])
  tree_copy = add_edge_into_index(tree_copy,edge_to_graft_on.dest,deleted_vertex,1,deleted_edge_ids[1])
  tree_copy = add_edge_into_index(tree_copy,grafted_subtree_root,deleted_vertex,1, graft_edge_id)

  return tree_copy

"""
  vertex_to_cut is an end point of cut_edge (cut_edge is not adjacent to root)
  grafting edge must be in component created after deleting cut_edge that does not contain the cut vertex
"""
def spr(tree,cut_edge_id,vertex_to_cut,grafting_edge):
  tree = deepcopy(tree)
  pruned_subtree,other_endpoint, deleted_edges = prune_subtree(tree,cut_edge_id,vertex_to_cut)
  grafted_tree = graft_subtree(pruned_subtree,grafting_edge,vertex_to_cut, other_endpoint, deleted_edges)
  return grafted_tree

def do_random_spr(tree):
  root_adjacent = []
  for i in range(len(tree.adjList[0])):
    e = tree.adjList[0,i]
    if e is None or e.id==-1:
      break
    root_adjacent.append(e.id)

  possible_cut_edges = []
  for i in range(len(tree.edges)):
    if i not in root_adjacent:
      possible_cut_edges.append(i)

  cut_edge = random.sample(possible_cut_edges,1)[0]

  def bfs(tree, root):
    parent = dict()
    q = collections.deque()
    q.append(root)
    parent[root] = -1

    edges = []

    while len(q)>0:
      u = q.popleft()
      for i in range(len(tree.adjList[u])):
        e = tree.adjList[u,i]
        if e is None or e.id==-1:
          break
        v = e.dest
        if v==u:
          v = e.src
        if v not in parent:
          parent[v] = u
          q.append(v)
          edges.append(e.id)
    return parent, edges

  parent, _ = bfs(tree, 0)
  e = tree.edges[cut_edge]
  vertex_to_cut = e.dest
  print(parent)
  if parent[vertex_to_cut] != e.src:
    vertex_to_cut = e.src

  ctree = deepcopy(tree)
  print(cut_edge, vertex_to_cut, tree.edges[cut_edge])
  pruned_subtree,other_endpoint, deleted_edges = prune_subtree(ctree,cut_edge,vertex_to_cut)
  print(pruned_subtree.edges)

  _, edges = bfs(pruned_subtree, 0)
  graft_edge_id = random.sample(edges,1)[0]
  print(graft_edge_id, pruned_subtree.edges[graft_edge_id])

  grafted_tree = graft_subtree(pruned_subtree,graft_edge_id,vertex_to_cut, other_endpoint, deleted_edges)
  return grafted_tree

def check_tree(tree_original, tree):
  assert(len(tree.edges) == tree.number_of_species*2-2)
  for i in range(len(tree.edges)):
    e = tree.edges[i]
    assert(e.id == i)
    u = e.src
    v = e.dest

    found_in_u = None
    for eu in tree.adjList[u]:
      if eu:
        if eu.id == i:
          if found_in_u == None:
            found_in_u = eu
          else:
            assert(False)
    assert(found_in_u != None)
    assert(found_in_u.src == u or found_in_u.dest == u)
    assert(found_in_u.src == v or found_in_u.dest == v)


    found_in_v = None
    for ev in tree.adjList[v]:
      if ev:
        if ev.id == i:
          if found_in_v == None:
            found_in_v = ev
          else:
            assert(False)
    assert(found_in_v != None)
    assert(found_in_v.src == u or found_in_v.dest == u)
    assert(found_in_v.src == v or found_in_v.dest == v)

  assert(len(tree.adjList) == tree.number_of_species*2-1)
  assert(len(tree.adjList) == tree.n)

  for i in range(tree.n):
    assert(len(tree.adjList[i]) == 4)
    cnt = 0
    for e in tree.adjList[i]:
      if e:
        assert(0 <= e.id < len(tree.edges))
        assert(0 <= e.src <tree.n)
        assert(e.src == i or e.dest == i)
        assert(e.src != e.dest)

        cnt+=1
    assert(cnt == tree.degree[i])
    assert(tree.adjList[i,-1] == None)

    flag = False
    for e in tree.adjList[i]:
      if e:
        assert(flag == False)
      else:
        flag = True
  assert(tree.degree[0] == 2)
  for i in range(1,tree.number_of_species+1):
    assert(tree.degree[i] == 1)
  for i in range(tree.number_of_species+1,tree.n):
    assert(tree.degree[i] == 3)
  for i in range(1,tree.number_of_species+1):
    assert(np.sum(tree.kmer_count[i] - tree_original.kmer_count[i]) == 0)

def generate_trees(original_tree, iteration, max_sprs):
  trees = [ [None for j in range(100)] for i in range(10) ]
  for i in range(max_sprs):
    for j in range(iteration):
      trees[i][j] = deepcopy(original_tree)
      for _ in range(i+1):
        trees[i][j] = do_random_spr(trees[i][j])
        check_tree(original_tree, trees[i][j])

  return trees

def get_rf_dist_dict(tree_original, trees, max_rf_distance, max_cnt):
  ans = { i : list() for i in range(max_rf_distance + 1) }
  for row in trees:
    for random_tree in row:
      dist = pt.rf_distance(tree_original, random_tree)
      if dist <= max_rf_distance and len(ans.get(dist)) < max_cnt:
        ans[dist].append(random_tree)

  return ans

def save_trees(outfile_prefix, tree_dict):
  for dist, trees in tree_dict.items():
    for i in range(len(trees)):
      pt.dump_to_file(trees[i], f'{outfile_prefix}{dist}-{i}-tree.txt')
      print(f'{outfile_prefix}{dist}-{i}-tree.txt')

def main():
  files = os.listdir(INPUT_TOPOLOGY_DIRECTORY)
  for file in files:
    original_treefile = os.path.join(INPUT_TOPOLOGY_DIRECTORY, file)
    original_tree = pt.load_from_file(original_treefile, only_leaves=False)
    trees = generate_trees(original_tree, ITERATIOINS, MAX_SPRS)

    tree_dicts = get_rf_dist_dict(original_tree, trees, MAX_RF_DISTANCE, MAX_CNT)

    outfile_prefix = os.path.join(OUTPUT_DIRECTORY, f'{file}-')
    save_trees(outfile_prefix, tree_dicts)

if __name__ == '__main__':
  main()
