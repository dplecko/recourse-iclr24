
from causallearn.search.ConstraintBased.FCI import fci
import numpy as np

heloc = r.data_top

G, edges = fci(heloc.to_numpy())

edge_list = []
for edge in edges:
    # only orient definitely direct edges
    if 1 in [x.value for x in edge.properties] or 2 in [x.value for x in edge.properties]:
        left_node = int(edge.get_node1().name[1:]) - 1
        right_node = int(edge.get_node2().name[1:]) - 1
        edge_list.append((left_node, right_node))

G_adj = {}
G_adj_par = {}
for i in range(len(heloc.columns)):
    G_adj[i] = []
    G_adj_par[i] = []

for (i, j) in edge_list:
    G_adj[i].append(j)
    G_adj_par[j].append(i)

G = np.zeros((len(heloc.columns), len(heloc.columns)))
for (i, j) in edge_list:
    G[i][j] = 1
