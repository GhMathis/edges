# -*- coding: utf-8 -*-
"""
Created on Tue May  9 17:52:49 2023

@author: Mathis
"""

from manim import *
from manim.utils.color import Colors
import numpy as np
import networkx as nx

class Networks(Scene):
    
    def adj_to_list(self,adj_arr):
        """Converte an adjacency matrix in the form of list of vertices and edges
        vertices = list(v1, v2, v3,...)
        edges = list((v1, v2), (v3, v1),...) """
        
        ##### vertices conversion
        self.vertices = np.arange(adj_arr.shape[0])
        ID = np.triu_indices(adj_arr.shape[0])
        
        ##### Edges conversion
        adj_arr[ID] = 0
        self.edges_temp = np.where(adj_arr == 1)
        self.edges = [(e1,e2) for e1,e2 in zip(self.edges_temp[0],self.edges_temp[1])]
    
    def egde_color(self,decompo_arr):
        ID = np.triu_indices(decompo_arr.shape[0],1)
        decompo_arr[ID] = 0
        pos_edges =  np.where(decompo_arr > 0)
        neg_edges =  np.where(decompo_arr < 0)
        egde_col_pos = {(e1,e2):{"stroke_color": Colors.blue.value} for e1,e2 in zip(pos_edges[0], pos_edges[1])}
        egde_col_pos2 = {(e1,e2):{"stroke_color": Colors.red.value} for e1,e2 in zip(neg_edges[0], neg_edges[1])}
        egde_col_pos = egde_col_pos| egde_col_pos2
        
        return(egde_col_pos)
    test = nx.gnp_random_graph(10,0.1)
    test.nodes
    test.edges
    def construct(self):
        adj_arr = np.array([[0, 1, 0, 1, 0, 1, 0, 1, 1, 0],
                            [1, 0, 0, 1, 0, 1, 1, 0, 1, 1],
                            [0, 0, 0, 0, 0, 0, 1, 1, 0, 0],
                            [1, 1, 0, 0, 0, 0, 1, 1, 1, 1],
                            [0, 0, 0, 0, 0, 1, 1, 0, 1, 1],
                            [1, 1, 0, 0, 1, 0, 1, 0, 0, 0],
                            [0, 1, 1, 1, 1, 1, 0, 0, 1, 1],
                            [1, 0, 1, 1, 0, 0, 0, 0, 0, 0],
                            [1, 1, 0, 1, 1, 0, 1, 0, 0, 1],
                            [0, 1, 0, 1, 1, 0, 1, 0, 1, 0]])
        
        spectra = np.linalg.eig(adj_arr)
        ##### Setup the data for visualisations
        self.adj_to_list(adj_arr)
        
        
        G_decomposition = lambda val, vec: np.dot(vec.T, vec)*np.exp(val)
        list_g_decomp = [G_decomposition(val, vec[np.newaxis]) for val, vec
                 in zip(spectra[0], spectra[1])]
        list_g_decomp[0] == list_g_decomp[0].T

        edge_col = self.egde_color(list_g_decomp[0])
        
        ##### Object creation
        #lyt = nx.bipartite_layout(self.vertices, self.vertices[:2])

        meta_ntw = [Graph(vertices = self.vertices,
                         edges = self.edges,
                         labels =True,
                         layout = "kamada_kawai",
                         edge_config = self.egde_color(g_decomp)).scale(0.5) for g_decomp
                    in list_g_decomp[:5]]
        meta_ntw2 = [Graph(vertices = self.vertices,
                         edges = self.edges,
                         labels =True,
                         layout = "spectral",
                         edge_config = self.egde_color(g_decomp)).scale(0.5) for g_decomp
                    in list_g_decomp[:5]]

        r1 = VGroup(*meta_ntw).arrange()
        r2 = VGroup(*meta_ntw2).arrange()
        
        self.play(LaggedStart(Create(VGroup(r1,r2).arrange(direction=DOWN)),
                              run_time=3))
        self.wait(2)
    
        
#python -m manim -pql graph_anime.py Networks