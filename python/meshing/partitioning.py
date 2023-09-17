import numpy as np
import pymetis, igl
from sklearn.cluster import (spectral_clustering, 
                            MeanShift, 
                            AffinityPropagation,
                            estimate_bandwidth)
                            


def spectralClustering(sf, n_parts):
    adjacency_matrix = igl.adjacency_matrix(sf)
    labels = spectral_clustering(
        adjacency_matrix, 
        n_components=n_parts)
    return labels

def  affinityPropagation(sf, n_parts):
    adjacency_matrix = igl.adjacency_matrix(sf)
    af = AffinityPropagation().fit(adjacency_matrix)
    cluster_centers_indices = af.cluster_centers_indices_
    labels = af.labels_
    return labels

def  meanShiftPartitioning(sf, n_parts):
    adjacency_matrix = igl.adjacency_matrix(sf)
    bandwidth = estimate_bandwidth(adjacency_matrix.toarray(), quantile=0.185, n_samples=600)
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=False)
    ms.fit(adjacency_matrix.toarray())
    labels = ms.labels_
    return labels

def metisPartitioning(sf, n_parts):
    adjacency_list = igl.adjacency_list(sf)
    n_cuts, membership = pymetis.part_graph(
        n_parts, 
        adjacency=adjacency_list)
    return membership

