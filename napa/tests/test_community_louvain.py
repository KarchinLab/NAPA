import networkx as nx
import napa.net.community_louvain as ml
from itertools import islice

g = nx.barbell_graph(4,5)
h = nx.barbell_graph(6,5)

mapping = zip(range(2*6+5), range(2*4+5+1,(2*6+5)+(2*4+5+1)))
print mapping 
mapping = dict(mapping)
mapping[8] = 7

h = nx.relabel_nodes(h,mapping, copy=False)
f = nx.compose(g, h)
print f.edges()
print 
dendo_dict = {r:ml.generate_dendrogram(graph = f, 
                                       weight = 'weight', resolution = r) \
              for r in [0.5, 1., 2., 4., 100., 1000]}
for r in sorted(dendo_dict.keys()):
    print r, dendo_dict[r]
    # print ml.partition_at_level(dendo_dict[r], level = len(dendo_dict[r])-1)
    best_part = ml.best_partition(f, resolution = r)
    print ml.modularity(best_part, f), best_part 
    print
                    

