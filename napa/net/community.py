# -*- coding: utf-8 -*-
'''
License
=======
NetworkX is distributed with the BSD license.

::

   Copyright (C) 2004-2016, NetworkX Developers
   Aric Hagberg <hagberg@lanl.gov>
   Dan Schult <dschult@colgate.edu>
   Pieter Swart <swart@lanl.gov>
   All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright 
  notice, this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright 
notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

* Neither the name of the NetworkX Developers nor the names of its
contributors may be used to endorse or promote products derived from 
this software without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
\"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=======
'''

'''
Copied from networkx 2.0 development version:
http://pelegm-networkx.readthedocs.io/en/latest/_modules/networkx/\
algorithms/community/centrality.html#girvan_newman
#-*- coding: utf-8 -*-  line needs to be first line
'''


import networkx as nx

__all__ = ['girvan_newman', 'asyn_lpa_communities']

# -*- coding: utf-8 -*-
# centrality.py - functions for computing communities using centrality notions
#
# Copyright 2015, 2016 NetworkX developers.
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; 
# see LICENSE.txt for more information.
"""Functions for computing communities based on centrality notions."""

def most_central_edge(G):
    centrality = nx.edge_betweenness_centrality(G, 
                                                weight='weight')
    return max(centrality, key=centrality.get)

def girvan_newman(G, most_valuable_edge=most_central_edge):
    """Finds communities in a graph using the Girvan–Newman method.

    Parameters
    ----------
    G : NetworkX graph

    most_valuable_edge : function
        Function that takes a graph as input and outputs an edge. The
        edge returned by this function will be recomputed and 
        removed at each iteration of the algorithm.

        If not specified, the edge with the highest
        :func:`networkx.edge_betweenness_centrality` will be used.

    Returns
    -------
    iterator
        Iterator over tuples of sets of nodes in `G`. 
        Each set of node is a community, 
        each tuple is a sequence of communities at a
        particular level of the algorithm.

    Examples
    --------
    To get the first pair of communities::

        >>> G = nx.path_graph(10)
        >>> comp = girvan_newman(G)
        >>> tuple(sorted(c) for c in next(comp))
        ([0, 1, 2, 3, 4], [5, 6, 7, 8, 9])

    To get only the first *k* tuples of communities, use
    :func:`itertools.islice`::

        >>> import itertools
        >>> G = nx.path_graph(8)
        >>> k = 2
        >>> comp = girvan_newman(G)
        >>> for communities in itertools.islice(comp, k):
        ...     print(tuple(sorted(c) for c in communities)) 
                #   doctest: +SKIP
        ...
        ([0, 1, 2, 3], [4, 5, 6, 7])
        ([0, 1], [2, 3], [4, 5, 6, 7])

    To stop getting tuples of communities once the number of   
    communities is greater than *k*, use 
    :func:`itertools.takewhile`::

        >>> import itertools
        >>> G = nx.path_graph(8)
        >>> k = 4
        >>> comp = girvan_newman(G)
        >>> limited = itertools.takewhile(lambda c: 
           len(c) <= k, comp)
        >>> for communities in limited:
        ...     print(tuple(sorted(c) for c in communities)) 
    # doctest: +SKIP
        ...
        ([0, 1, 2, 3], [4, 5, 6, 7])
        ([0, 1], [2, 3], [4, 5, 6, 7])
        ([0, 1], [2, 3], [4, 5], [6, 7])

    To just choose an edge to remove based on the weight::

        >>> from operator import itemgetter
        >>> G = nx.path_graph(10)
        >>> edges = G.edges()
        >>> nx.set_edge_attributes(G, 'weight', 
        ... {(u, v): v for u, v in edges})
        >>> def heaviest(G):
        ...     u, v, w = max(G.edges(data='weight'),      
        ...     key=itemgetter(2))
        ...     return (u, v)
        ...
        >>> comp = girvan_newman(G, most_valuable_edge=heaviest)
        >>> tuple(sorted(c) for c in next(comp))
        ([0, 1, 2, 3, 4, 5, 6, 7, 8], [9])

    To utilize edge weights when choosing an edge with, 
    for example, the highest betweenness centrality::

        >>> from networkx import edge_betweenness_centrality \
        ... as betweenness
        >>> def most_central_edge(G):
        ...     centrality = betweenness(G, weight='weight')
        ...     return max(centrality, key=centrality.get)
        ...
        >>> G = nx.path_graph(10)
        >>> comp = girvan_newman(G, 
            most_valuable_edge=most_central_edge)
        >>> tuple(sorted(c) for c in next(comp))
        ([0, 1, 2, 3, 4], [5, 6, 7, 8, 9])

    To specify a different ranking algorithm for edges, use the
    `most_valuable_edge` keyword argument::

        >>> from networkx import edge_betweenness_centrality
        >>> from random import random
        >>> def most_central_edge(G):
        ...     centrality = edge_betweenness_centrality(G)
        ...     max_cent = max(centrality.values())
        ...     # Scale the centrality values so they are 
        ...     # between 0 and 1,
        ...     # and add some random noise.
        ...     centrality = \
        ...     {e: c / max_cent for e, c in centrality.items()}
        ...     # Add some random noise.
        ...     centrality = {e: c + random() for e, c in       
        ...     centrality.items()}
        ...     return max(centrality, key=centrality.get)
        ...
        >>> G = nx.path_graph(10)
        >>> comp = girvan_newman(G, 
        most_valuable_edge=most_central_edge)

    Notes
    -----
    The Girvan–Newman algorithm detects communities by progressively
    removing edges from the original graph. The algorithm removes the
    "most valuable" edge, traditionally the edge with the highest
    betweenness centrality, at each step. 
    As the graph breaks down into pieces, the tightly knit community 
    structure is exposed and the result can be depicted as a 
    dendrogram.

    """
    # If the graph is already empty, simply return its connected
    # components.
    if G.number_of_edges() == 0:
        yield tuple(nx.connected_components(G))
        return
    # If no function is provided for computing the 
    # most valuable edge,
    # use the edge betweenness centrality.
    if most_valuable_edge is None:
        def most_valuable_edge(G):
            """Returns the edge with the 
            highest betweenness centrality
            in the graph `G`.

            """
            # We have guaranteed that the graph is non-empty, so this
            # dictionary will never be empty.
            betweenness = nx.edge_betweenness_centrality(G)
            return max(betweenness, key=betweenness.get)
    
    # The copy of G here must include the edge weight data.
    g = G.copy().to_undirected()
    
    # Self-loops must be removed because their removal 
    # has no effect on
    # the connected components of the graph.
    g.remove_edges_from(g.selfloop_edges())
    while g.number_of_edges() > 0:
        yield _without_most_central_edges(g, most_valuable_edge)



def _without_most_central_edges(G, most_valuable_edge):
    """Returns the connected components of the graph that 
    results from repeatedly removing the most "valuable" 
    edge in the graph.

    `G` must be a non-empty graph. 
    This function modifies the graph `G`
    in-place; that is, it removes edges on the graph `G`.

    `most_valuable_edge` is a function that takes 
    the graph `G` as input
    (or a subgraph with one or more edges of `G` removed) and 
    returns an edge. 
    That edge will be removed and this process will be repeated
    until the number of connected components in the graph increases.

    """
    original_num_components = nx.number_connected_components(G)
    num_new_components = original_num_components
    while num_new_components <= original_num_components:
        edge = most_valuable_edge(G)
        G.remove_edge(*edge)
        new_components = tuple(nx.connected_components(G))
        num_new_components = len(new_components)
    return new_components

'''
def girvan_newman(G, weight=None):
    """
    Find communities in graph using Girvan–Newman method.
    
    Parameters
    ----------
    G : NetworkX graph

    weight : string, optional (default=None)
       Edge data key corresponding to the edge weight.

    Returns
    -------
    List of tuples which contains the clusters of nodes.

    Examples
    --------
    >>> G = nx.path_graph(10)
    >>> comp = girvan_newman(G)
    >>> comp[0]
    ([0, 1, 2, 3, 4], [8, 9, 5, 6, 7])

    Notes
    -----
    The Girvan–Newman algorithm detects communities by progressively
    removing edges from the original graph. Algorithm removes edge 
    with the highest betweenness centrality at each step. 
    
    As the graph breaks down into pieces, the tightly
    knit community structure is exposed and result can be depicted 
    as a dendrogram.
    
    """
    g = G.copy().to_undirected()
    
    components = []
    while g.number_of_edges() > 0:
        _remove_max_edge(g, weight)
        components.append(tuple(list(H)
            for H in \
                nx.connected_component_subgraphs(g)))

    return components

def _remove_max_edge(G, weight=None):
    """
    Removes edge with the highest value on 
    betweenness centrality.

    Repeat this step until more connected components 
    than the connected components of the original 
    graph are detected.

    It is part of Girvan–Newman algorithm.

    :param G: NetworkX graph
    :param weight: string, optional (default=None) 
                   Edge data key corresponding
                   to the edge weight.
    """
    
    number_components = \
    nx.number_connected_components(G)
    
    while nx.number_connected_components(G) <= \
          number_components:
        betweenness = \
        nx.edge_betweenness_centrality(G, 
                                       weight=weight)
        
        max_value = max(betweenness.values())
        
        # Use a list of edges because G is 
        # changed in the loop
        for edge in list(G.edges()):
            if betweenness[edge] == max_value:
                G.remove_edge(*edge)

'''

#=========================================================#
# -*- coding: utf-8 -*-
#    Copyright (C) 2015
#    All rights reserved.
#    BSD license.
"""
Asynchronous label propagation algorithms for 
community detection.
"""

from collections import Counter
from collections import defaultdict
import random

#from networkx.utils import groups

def groups(many_to_one):
    """
    Converts a many-to-one mapping into a one-to-many mapping.
    `many_to_one` must be a dictionary whose keys and values are all
    :term:`hashable`.
    The return value is a dictionary mapping values from 
    `many_to_one` to sets of keys from `many_to_one` that have that 
    value.
    
    For example::
        >>> from networkx.utils import groups
        >>> many_to_one = {'a': 1, 'b': 1, 
                           'c': 2, 'd': 3, 'e': 3}
        >>> groups(many_to_one)  # doctest: +SKIP
        {1: {'a', 'b'}, 2: {'c'}, 3: {'d', 'e'}}
    
    """
    one_to_many = defaultdict(set)
    for v, k in many_to_one.items():
        one_to_many[k].add(v)
    return dict(one_to_many)


def asyn_lpa_communities(G, weight=None):
    """
    Returns communities in `G` as detected by asynchronous label
    propagation.

    The asynchronous label propagation algorithm is described in 
    [1]_. 
    The algorithm is probabilistic and the found communities may
    vary on different executions.

    The algorithm proceeds as follows. 

    After initializing each node with a unique label, the algorithm 
    repeatedly sets the label of a node to be the label that appears 
    most frequently among that node`s neighbors. 

    The algorithm halts when each node has the label that appears 
    most frequently among its neighbors. 
    
    The algorithm is asynchronous because each node is updated 
    without waiting for updates on the remaining nodes.

    This generalized version of the algorithm in [1]_ accepts edge 
    weights.

    Parameters
    ----------
    G : Graph

    weight : string
        The edge attribute representing the weight of an edge.
        If None, each edge is assumed to have weight one. 
        In this algorithm, the weight of an edge is used in 
        determining the frequency with which a label appears 
        among the neighbors of a node: a higher weight means the 
        label appears more often.

    Returns
    -------
    communities : dictionary
        Dictionary of communities given as sets of nodes.

    Notes
    ------
    Edge weight attributes must be numerical.


    References
    ----------
    .. [1] Raghavan, Usha Nandini, Réka Albert, 
           and Soundar Kumara. "Near linear time 
           algorithm to detect community structures 
           in large-scale networks." 
           Physical Review E 76.3 (2007): 036106.
    """

    labels = {n: i for i, n in enumerate(G)}
    cont = True
    while cont:
        cont = False
        nodes = list(G)
        random.shuffle(nodes)
        # Calculate the label for each node
        for node in nodes:
            if len(G[node]) < 1:
                continue

            # Get label frequencies. Depending on 
            # the order they are processed in 
            # some nodes with be in t and others in t-1, 
            # making the algorithm asynchronous.
            label_freq = \
            Counter({labels[v]: G.edge[v][node][weight]
                     if weight else 1 for v in G[node]})

            # Choose the label with the highest frecuency. 
            # If more than 1 label has the highest 
            # frecuency choose one randomly.
            max_freq = max(label_freq.values())

            best_labels = [label for label, freq in \
                           label_freq.items() \
                           if freq == max_freq]

            new_label = random.choice(best_labels)

            labels[node] = new_label

            # Continue until all nodes have a 
            # label that is better than other
            # neighbour labels (only one label 
            # has max_freq for each node).
            cont = cont or len(best_labels) > 1

    # TODO In Python 3.3 or later, 
    # this should be `yield from ...`.

    # return iter(groups(labels).values())
    return groups(labels).values()
