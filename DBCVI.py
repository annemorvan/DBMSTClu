# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 10:55:08 2017

@author: amorvan
"""

from __future__ import division

import igraph as ig
import sys
import numpy as np
from collections import deque



#def get_edges_list(G, weight = False):
#    """
#    Returns the list of edges of G in the form [ ('node1', 'node2'), ...].
#    If weight = True, add also the weight: [ ('node1', 'node2', w), ...]. 
#    """
#    if not weight:
#        return [ (G.vs[edge.source]["name"], G.vs[edge.target]["name"]) for edge in G.es]
#    else:
#        return [ (G.vs[edge.source]["name"], G.vs[edge.target]["name"], edge['weight']) for edge in G.es]
#
#
#
#def get_nodes_list(G):
#    """
#    Returns the list of nodes in the form [ 'node1', 'node2', ...]
#    """
#    return [node["name"] for node in G.vs] 



def evaluateCuts(T0, N, nbCluster, cluster_info, verbose = False):
    """
    Update size_source/target, disp_source/target and sep_source/target
    for each possible cut in the corresponding connected component in the mst.
    Update sep_node and ix_cluster for each node.
    """
    
    # each cluster is represented by one of its edge + its number of points + VC
    #T0D0: 28/02/2018
    node_ix, other_node_ix, old_size, old_VC, _ = cluster_info 
    if other_node_ix == None: # the cluster is a singleton, no cuts to evaluate!
        assert old_size == 1
#        T0.vs[node_ix]['ix_cluster'] = nbCluster
        if verbose:
            print('THE CONSIDERED CLUSTER IS A SINGLETON')
        return None

    eid = T0.get_eid(node_ix, other_node_ix)
    assert eid == T0.get_eid(other_node_ix, node_ix)

    marked = False
    
    # variables to return only the optimal cut
    current_optimal_diff_DBCVI = -2.0
    optimal_cut = eid
    
    l_edges_to_evaluate = deque()
    
    if verbose:
        print('')
        print('current processed edge (eid, node_ix, marked)', eid, node_ix, marked)
            
    source_eid_list = [incident_ix for incident_ix in T0.incident(node_ix) if (incident_ix != eid) ]
    target_eid_list = [incident_ix for incident_ix in T0.incident(other_node_ix) if (incident_ix != eid) ]
    
    # put the edge itself marked, and see the other side
    l_edges_to_evaluate.append( (eid, other_node_ix, True))
#    if verbose:
#        print('l_edges_to_evaluate after adding the edge itself marked with other side', l_edges_to_evaluate)
    
    # put incident edges from the other side
    for target_eid in target_eid_list:
        if T0.es[target_eid].source == other_node_ix:
            l_edges_to_evaluate.append( (target_eid, T0.es[target_eid].target, False))    
        else:
            l_edges_to_evaluate.append( (target_eid, T0.es[target_eid].source, False))   
#    if verbose:
#        print('l_edges_to_evaluate after adding the incident edges from the other side', l_edges_to_evaluate)
    
    # put the edge itself marked, for the initial side
    l_edges_to_evaluate.append( (eid, node_ix, True)) 
#    if verbose:
#        print('l_edges_to_evaluate after the edge itself marked with the initial side', l_edges_to_evaluate)
    
    # put incident edges from the initial side
    for source_eid in source_eid_list:
        if T0.es[source_eid].source == node_ix:
            l_edges_to_evaluate.append( (source_eid, T0.es[source_eid].target, False))    
        else:
            l_edges_to_evaluate.append( (source_eid, T0.es[source_eid].source, False)) 
    if verbose:
        print('l_edges_to_evaluate after adding the incident edges from the initial side', l_edges_to_evaluate)            
    
    while len(l_edges_to_evaluate) != 0 :
        
        eid, node_ix, marked = l_edges_to_evaluate.pop()
        
        if verbose:
            print('')
            print('current processed edge (eid, node_ix, marked)', eid, node_ix, marked)
        
        # if the edge has not already been seen
        if not marked:  
            
            if verbose:
                print('the edge was not marked')
            
            incident_eid_list = [incident_ix for incident_ix in T0.incident(node_ix) if (incident_ix != eid) ]
            
            # put the edge itself marked, for the other side (not in priority)
            if T0.es[eid].source == node_ix:
                other_node_ix = T0.es[eid].target
            else:
                other_node_ix = T0.es[eid].source
                
            #l_edges_to_evaluate.insert(0,  (eid, other_node_ix, True))
            l_edges_to_evaluate.appendleft( (eid, other_node_ix, True) )
            
            if verbose:
                print('l_edges_to_evaluate after adding the edge itself at the beginning, other side', l_edges_to_evaluate)  

            # put the edge itself marked, initial side
            l_edges_to_evaluate.append( (eid, node_ix, True))  
            if verbose:
                print('l_edges_to_evaluate after adding the edge itself initial side', l_edges_to_evaluate)   
                
            # put incident edges 
            for incident_eid in incident_eid_list:
                if T0.es[incident_eid].source == node_ix:
                    l_edges_to_evaluate.append( (incident_eid, T0.es[incident_eid].target, False))    
                else:
                    l_edges_to_evaluate.append( (incident_eid, T0.es[incident_eid].source, False)) 

            if verbose:
                print('l_edges_to_evaluate', l_edges_to_evaluate)            
    
        else:
            if verbose:
                print('the edge is marked')
            
            # affectation to a cluster
            T0.vs[node_ix]['ix_cluster'] = nbCluster
            
            if node_ix == T0.es[eid].source:
                side = 'source'
            else:
                side = 'target'
            
            incident_eid_list = [incident_ix for incident_ix in T0.incident(node_ix) if (incident_ix != eid) ]
                         
            if verbose:
                print('incident_eid_list', incident_eid_list)
            
            # for updating disp
            incident_disp_list = []
            incident_weights_list = []

            # for updating sep
            w = T0.es[eid]['weight']
            incident_sep_list = []
            
            # for updating size
            T0.es[eid]['size_' + side] = 0
            
            for incident_eid in incident_eid_list:
                
                if T0.es[incident_eid].source == node_ix:
                    incident_side = 'target'
                else:
                    incident_side = 'source'
                
                # update size
                T0.es[eid]['size_' + side] += T0.es[incident_eid]['size_' + incident_side]
                
                # for updating disp
                incident_disp_list.append(T0.es[incident_eid]['disp_' + incident_side])
                incident_weights_list.append(T0.es[incident_eid]['weight'])
                
                # for updating sep
                incident_sep_list.append(T0.es[incident_eid]['sep_node_' + incident_side])
            
            # update size
            T0.es[eid]['size_' + side] += 1
                    
            # update disp
            if verbose:
                print('incident_disp_list', incident_disp_list)
                print('incident_weights_list', incident_weights_list)
            
            if len(incident_disp_list) != 0:
                
                # update disp
                maxi_disp = max(incident_disp_list)
                maxi_weights = max(incident_weights_list)
                T0.es[eid]['disp_' + side] = max(maxi_disp, maxi_weights)    
            
                # update sep
                mini_sep = min(incident_sep_list)
                T0.es[eid]['sep_node_' + side] = min(mini_sep, T0.vs[node_ix]['sep_node'])
                
            else:
                # update sep
                T0.es[eid]['sep_node_' + side] = T0.vs[node_ix]['sep_node']
                T0.es[eid]['disp_' + side] = 0.0 
                
            
            # update VC
            true_sep = min(T0.es[eid]['sep_node_' + side], w)
#            T0.es[eid]['sep_' + side] = true_sep
            if verbose:
                print('true_sep', true_sep)
            T0.es[eid]['VC_' + side] = (true_sep - T0.es[eid]['disp_' + side]) / max (true_sep, T0.es[eid]['disp_' + side])
            assert T0.es[eid]['VC_' + side] >= -1.001
            assert T0.es[eid]['VC_' + side] <= 1.001
            if verbose:
                print('size_' + side, 'disp_' + side, 'sep_node_' + side, T0.es[eid]['size_' + side], T0.es[eid]['disp_' + side], T0.es[eid]['sep_node_' + side] )            
   
            # update side_ok
            T0.es[eid][ side + '_ok'] = True   
   
            # update diff_DBCVI
            new_partial_DBCVI = (T0.es[eid]['size_source'] * T0.es[eid]['VC_source'] + T0.es[eid]['size_target'] * T0.es[eid]['VC_target']) / N 
            if verbose:
                print('T0.es[eid][size_source]', T0.es[eid]['size_source'])
                print('T0.es[eid][VC_source]', T0.es[eid]['VC_source'])
                print('T0.es[eid][size_target]', T0.es[eid]['size_target'])
                print('T0.es[eid][VC_target]', T0.es[eid]['VC_target'])
                print('new_partial_DBCVI', new_partial_DBCVI)
                print('(old_size * old_VC)/N', (old_size * old_VC)/N)
            T0.es[eid]['diff_DBCVI'] = new_partial_DBCVI - (old_size * old_VC)/N
            if verbose:
                print('T0.es[eid][diff_DBCVI]', T0.es[eid]['diff_DBCVI'])

            assert T0.es[eid]['diff_DBCVI'] >= -2.001
            assert T0.es[eid]['diff_DBCVI'] <= 2.001
            
            if T0.es[eid]['source_ok'] and T0.es[eid]['target_ok'] and T0.es[eid]['diff_DBCVI'] > current_optimal_diff_DBCVI:
                optimal_cut = eid
                current_optimal_diff_DBCVI = T0.es[eid]['diff_DBCVI']
                    
    # return the optimal cut for the considered connected component
    return optimal_cut

   
            
def DBMSTClu(T0, initial_DBCVI = -1.0, stop_criterion = 1.0, max_nb_clusters = sys.maxsize, verbose = False, path_save_cuts = None):
    """
    DBMSTClu algorithm.
    Parameters:
    T0: MST (igraph graph)
    Begin with the MST and find at each iteration the best cut
    regarding the DBCVI.
    
    With igraph, to obtain:
    * the list of edges                                     : T0.es
    * the list of nodes                                     : T0.vs
    * for each node, the list of incident edges to this node: T0.incident(node) or T0.incident('node')
    * for each node, the list of its direct neighbors       : T0.neighbors(node) or T0.incident('node')
    
    Lots of information to store about nodes and edges are attributes of nodes or edges because of the continuous indexing
    of nodes and edges.
    
    THE CASE WITH DISCONNECTED MST HAS NOT BEEN TESTED YET.
    """    
    assert T0.vs['name']
    assert T0.es['weight']
    assert all(w > 0.0 for w in T0.es['weight'])
    #TODO:
    assert T0.is_connected() #THE CASE WITH DISCONNECTED MST HAS NOT BEEN TESTED YET.
    
    # initialization
    N = T0.vcount() # number of nodes
    M = T0.ecount() # number of edges
    if T0.is_connected():
        assert M == (N-1)
    assert M >= 1

    # initialization of the associate index of each edge to clusters
    T0.vs['ix_cluster'] = [0]*N # to modify wrt to the previous stage

    # initialization of other attributes of edges and nodes
    T0.vs['sep_node'] = [sys.maxsize]*N  
    
    T0.es['sep_node_source'] = [sys.maxsize]*M
    T0.es['disp_source'] = [0.0]*M
    T0.es['size_source'] = [0]*M
    T0.es['VC_source'] = [-1.0]*M
    T0.es['source_ok'] = [False]*M

    T0.es['sep_node_target'] = [sys.maxsize]*M
    T0.es['disp_target'] = [0.0]*M
    T0.es['size_target'] = [0]*M
    T0.es['VC_target'] = [-1.0]*M
    T0.es['target_ok'] = [False]*M
    
    # diff_DBCVI
    T0.es['diff_DBCVI'] = [0.0]*M

    
    # initialization of the list of clusters 
    # the list has the size of the number of clusters
    # each cluster is represented by one of its edge (node1, node2) + its number of points + VC + its indexing in the cluster_infos
    # the edge can not be represented as eid because only nodes are not deleted (continuous indexing for nodes and edges with igraph!)
    node_ix, other_node_ix = T0.es[0].source, T0.es[0].target
    clusters_info = [ (node_ix, other_node_ix, N, -1.0, 0) ] 

    # initialization of the number of cuts
    nbCluster = 0

    # initialization of the list of cuts
    cuts = []
    # initialization of the list of optimal cuts per cluster
    optimal_cuts_list = []    
    
    # initialization of the DBCVI
    final_DBCVI = initial_DBCVI
    
    if verbose:
        print('')
        print('stop criterion DBCVI', stop_criterion)
        print('')
    
    while len(cuts) < M : 
    #while final_DBCVI < stop_criterion and len(cuts) < (max_nb_clusters - 1): # while we need cuts
        
        # look for an edge getting cut in the last two formed clusters or in the initial big cluster 
        len_clusters_info = len(clusters_info)
        
        for i in range( max(len_clusters_info - 2, 0), len_clusters_info):
            if verbose or True:
                print('evaluation of ' + str(i) + 'th cluster: ' + str(clusters_info[i]))
            optimal_cut_cluster_eid = evaluateCuts(T0, N, nbCluster, clusters_info[i], verbose) # will update ix_cluster
            if optimal_cut_cluster_eid != None: # if the evaluted cluster is not a singleton...
                if verbose or True:
                    print('an optimal cut (eid) added for cluster ' + str(i) + ': ' + str(optimal_cut_cluster_eid))
                optimal_cuts_list.append( (T0.es[optimal_cut_cluster_eid].source, T0.es[optimal_cut_cluster_eid].target))            
            nbCluster += 1
        
        if verbose or True:
            print('before cut, optimal_cuts_list', optimal_cuts_list)

        if len(optimal_cuts_list) == 0.0 :
            raise ValueError("It would mean that all optimal cuts per cluster build singletons => impossible since len(cuts) < M.")
        if verbose or True:
            for node_eid_optimal, other_node_eid_optimal in optimal_cuts_list:
                print('considered cut', node_eid_optimal, other_node_eid_optimal)
                print('considered diff_DBCVI', T0.es[T0.get_eid(node_eid_optimal, other_node_eid_optimal)]['diff_DBCVI']) 

        node_eid_max, other_node_eid_max = optimal_cuts_list[np.argmax([ T0.es[ T0.get_eid(node_eid_optimal, other_node_eid_optimal)]['diff_DBCVI'] for (node_eid_optimal, other_node_eid_optimal) in optimal_cuts_list])]
        eid_max = T0.get_eid(node_eid_max, other_node_eid_max)
        assert eid_max == T0.get_eid(other_node_eid_max, node_eid_max)
        diff_DBCVI_max = T0.es[eid_max]['diff_DBCVI']
            
        if diff_DBCVI_max >= 0.0: # if the cut found improves the current DBCVI
            
            # perform the cut + update relevant structures not updated in evaluateCuts
            edge = T0.es[eid_max]
            node1_ix, node2_ix = edge.source, edge.target
            node1, node2 = T0.vs[node1_ix]['name'], T0.vs[node2_ix]['name'] 
            w = T0.es[eid_max]['weight'] 
            
            # update cuts
            cuts.append( (node1, node2, w) )
            #print('cuts', len(cuts), cuts[-1])
            
            # update final_DBCVI
            if verbose:
                print('')
                print('-------> diff_DBCVI_max, nbCluster', diff_DBCVI_max, nbCluster)
                print('associated cut (node1, node2, w)', node1, node2, w)           
                print('final DBCVI before', final_DBCVI)
                
            final_DBCVI += diff_DBCVI_max
            if verbose or True:
                print('final DBCVI after', final_DBCVI)
            assert final_DBCVI >= -1.001
            assert final_DBCVI <= 1.001
            
            # update separation for source and target points
            T0.vs[node1_ix]['sep_node'] = min(T0.vs[node1_ix]['sep_node'] , w)
            T0.vs[node2_ix]['sep_node'] = min(T0.vs[node2_ix]['sep_node'] , w)
            
            # get the info for clusters_info before cutting te edge
            size_source, VC_source = edge['size_source'], edge['VC_source']
            size_target, VC_target = edge['size_target'], edge['VC_target'] 
            if verbose:
                print('size_source, size_target', size_source, size_target)
            
            # update clusters_info = [ (node_ix, other_node_ix, size, partial_DBCVI, noCluster) ]
            cluster_id_to_find = T0.vs[node1_ix]['ix_cluster']
            if verbose:
                print('clusters_info', clusters_info)
                print('cluster_id_to_find', cluster_id_to_find)
            to_delete = [noCluster for _,_,_,_,noCluster in clusters_info].index( cluster_id_to_find )
            if verbose:
                print('to_delete', to_delete)
            del clusters_info[to_delete] 
                        
            # finally perform definitly the cut
            T0.delete_edges(eid_max) 
            del optimal_cuts_list[ optimal_cuts_list.index( (node_eid_max, other_node_eid_max))]
            if verbose or True:
                print('after cut: updated optimal_cuts_list', optimal_cuts_list)
            
            # /!\ To do after having cut the edge because otherwise the indexing changes!
            # take a source edge of the cut edge (incident to source), here the first one
            source_eid_list = T0.incident(node1_ix)
            if len(source_eid_list) != 0:
                source_eid = T0.incident(node1_ix)[0]
                clusters_info.append( (T0.es[source_eid].source, T0.es[source_eid].target, size_source, VC_source, nbCluster) )                
            else:
                # SINGLETON CLUSTER!!!
                clusters_info.append( (node1_ix, None, 1, 1.0, nbCluster) )   
                T0.vs[node1_ix]['ix_cluster'] = nbCluster
            # take a right edge of the cut edge (incident to target) 
            target_eid_list = T0.incident(node2_ix)
            if len(target_eid_list) != 0:
                target_eid = T0.incident(node2_ix)[0] 
                clusters_info.append( (T0.es[target_eid].source, T0.es[target_eid].target, size_target, VC_target, nbCluster + 1) )
            else:
                # SINGLETON CLUSTER!!!
                clusters_info.append( (node2_ix, None, 1, 1.0, nbCluster + 1) )
                T0.vs[node2_ix]['ix_cluster'] = nbCluster + 1
            # nbCluster and not T0.vs[node1_ix]['ix_cluster'] or T0.vs[node2_ix]['ix_cluster'] because it has not been updated yet.
            if verbose or True:
                print('after cut, updated clusters_info', clusters_info)
                print('\n')

            
            
        else:
            break
    
    if verbose:    
        print('cuts in DBCVI function')
        print(cuts)
        if len(cuts) == max_nb_clusters:
            print('*** MAXIMUM NUMBER OF CUTS REACHED ***', len(cuts))
        print('')
    return final_DBCVI, cuts, T0.vs['ix_cluster']
    

