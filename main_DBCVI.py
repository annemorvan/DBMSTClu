#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 12:56:34 2018

@author: amorvan
"""
from __future__ import division
import numpy as np
import os
from sklearn import datasets
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import euclidean_distances
from faker import Faker
import matplotlib.pyplot as plt
import igraph as ig

import DBCVI


if __name__ == '__main__':

    if not os.path.exists('./figures/moons/'):
        os.makedirs('./figures/moons/')     
    
    print_labels = True
    print_artificial_colors = True    
    # parameters:    
    seed = 0
    n_samples = 100
    assert n_samples % 2 == 0
    add_name = ''    
    

    X, y = datasets.make_moons(n_samples, shuffle = False, noise = 0.09, random_state = seed)
    assert X.shape == (n_samples, 2)
    
    plt.subplot(111)
    #plt.title(name, size=18)
    colors = ['b', 'g', 'r', 'c', 'm', 'k']
    plt.scatter(X[:, 0], X[:, 1], s = 8, color = [colors[i] for i in y])
    plt.xlim(-1.0, 1.0)
    plt.ylim(-1.25, 1.25)    
    plt.grid()
    plt.axis('equal')
#    plt.xticks(())
#    plt.yticks(())
    plt.show()    
    
    
    # generate the matrix of distances
    DistEuclid = euclidean_distances(X)
    assert DistEuclid.shape == (n_samples, n_samples)
    dist_max = np.max(DistEuclid)
    DistEuclid /= dist_max
    print('before')
    print(DistEuclid)
    
    A = DistEuclid    
    n_edges = int(np.count_nonzero(A) / 2)
    print('Number of edges', n_edges)
    ratio = n_edges / n_samples
    print('Ration between the number of edges and the number of nodes', ratio)

    # Create graph, A_triu.astype(bool).tolist() or (A_triu / A_triu).tolist() can also be used.
    colors = ['#0000FF', '#00FF00', '#FF0000', '#FFFF00', '#FF00FF', '#00FFFF', '#00FF00', '#008080', '#800000', '#808000', '#800080', '#808080']
    A_triu = np.triu(A, 0)
    print('A_triu')
    print(A_triu)
    weights_list = A_triu[A_triu.nonzero()]   
    assert len(weights_list) == n_edges
    g = ig.Graph.Adjacency((A_triu > 0.0).tolist(), mode = ig.ADJ_UNDIRECTED) 
    # Add edge weights and node labels.
    g.es['weight'] = weights_list
    g.vs['name'] = [str(i) for i in xrange(n_samples)]  
    if print_labels:
        g.vs['label'] = [str(i) for i in xrange(n_samples)]     
    g.vs['color'] = ['#0000FF']*(n_samples // 2) + ['#00FF00']*(n_samples // 2)
    assert len(g.es) == n_edges
    assert n_edges >= (n_samples - 1)

    # compute mst
    mst = g.spanning_tree( weights = g.es['weight'])
    print('len(mst.es)', len(mst.es))
    assert len(mst.es) == n_samples - 1
    
    # Plot the graph 
    visual_style = {}
    # Scale vertices based on degree
    visual_style["vertex_size"] = [20] * n_samples
    # Don't curve the edges
    visual_style["edge_curved"] = False    
    # Set bbox and margin
    visual_style["bbox"] = (800,800)
    visual_style["margin"] = 100 
    visual_style["layout"] = [ (X[i,0], X[i,1]) for i in xrange(n_samples)]
    ig.plot(mst, 'figures/moons/moons_' + add_name + 'mst_n' + str(n_samples) + '_m' + str(n_edges) + '.pdf', **visual_style)
    ig.plot(g, 'figures/moons/moons_' + add_name + 'graph_n' + str(n_samples) + '_m' + str(n_edges) + '.pdf', **visual_style)


    # APPLY KMEANS
    kmeans = KMeans(n_clusters = 2, random_state = 0).fit(X)  
    mst.vs['color'] = [ colors[i] for i in list(kmeans.labels_)]
    ig.plot(mst, 'figures/moons/moons_' + add_name + 'mst_n' + str(n_samples) + '_m' + str(n_edges) + '_Kmeans.pdf', **visual_style)

    
    # APPLY CLUSTERING DBMSTClu
    index_DBMSTClu, DBMSTClu_cuts, y_pred_ix = DBCVI.DBMSTClu(mst, verbose = False)
    print('final DBMSTClu index', index_DBMSTClu)
    print('final DBMSTClu cuts', DBMSTClu_cuts)
    print('final y_pred_ix', y_pred_ix)
    n1, n2 = DBMSTClu_cuts[0][0], DBMSTClu_cuts[0][1]
    mst.add_edges([(n1, n2)])
    eid = mst.get_eid(n1, n2)    
    unique, counts = np.unique( np.asarray(y_pred_ix), return_counts = True)
    unique, counts = list(unique), list(counts)
    # generate random colors
    if not print_artificial_colors:
        fake = Faker()
        colors = [fake.hex_color() for _ in xrange(len(unique))]
    
    unique_copy = list(unique)
    unique.sort( key = lambda x: counts[unique_copy.index(x)], reverse = True)
    mst.vs['color'] = [colors[ unique.index(i)] for i in y_pred_ix]
    ig.plot(mst, 'figures/moons/moons_' + add_name + 'mst_n' + str(n_samples) + '_m' + str(n_edges) + '_DBMSTClu.pdf', **visual_style)




