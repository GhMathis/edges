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
    def verte_color(self,eign_vec):

        pos_verte =  np.where(eign_vec > 0)
        neg_verte =  np.where(eign_vec < 0)
        
        verte_col_pos = {v:{"fill_color": Colors.blue.value} for v in pos_verte[0]}
        verte_col_pos2 = {v:{"fill_color": Colors.red.value} for v in neg_verte[0]}
        verte_col_pos = verte_col_pos| verte_col_pos2
        
        return(verte_col_pos)
    def construct(self):
        cbPalette = ["#999999", "#E69F00", "#56B4E9",
                      "#009E73", "#F0E442", "#0072B2",
                      "#D55E00", "#CC79A7"]
        self.camera.background_color = "#FFFFFF"
        
        adj_arr = np.array([[ 0, 0, 0, 1, 0],
                            [ 0, 0, 0, 1, 1],
                            [ 0, 0, 0, 1, 0],
                            [ 1, 1, 1, 0, 0],
                            [ 0, 1, 0, 0, 0]])
        
        
        
        
        ##### Math, matrix and spectre computation
        spectra = np.linalg.eig(adj_arr)
        self.adj_to_list(adj_arr)
        
        G_decomposition = lambda val, vec: np.dot(vec.T, vec)*np.exp(val)
        list_g_decomp = [G_decomposition(val, vec[np.newaxis]) for val, vec
                 in zip(spectra[0], spectra[1])]
        list_g_decomp[0] == list_g_decomp[0].T
        
        ##### aesthetics
        partitions = [[i for i in range (3)], [i for i in range (3,5)]]
        
        vertex_config = {v:{"fill_color": cbPalette[1]} for v in 
                          [i for i in range (5)]}
        

        ##### Object creation
        #lyt = nx.bipartite_layout(self.vertices, self.vertices[:2])
        
        meta_ntw = [Graph(vertices = self.vertices,
                         edges = self.edges,
                         labels =True,
                         layout = "partite",
                         partitions = partitions ,
                         vertex_config = vertex_config,#self.verte_color( vec*np.exp(val))
                         edge_config = self.egde_color(g_decomp)) for val, vec, g_decomp
                    in zip(spectra[0][:5], spectra[1][:5],list_g_decomp[:5])]


        r1 = VGroup(*meta_ntw)

        
        self.add(r1.arrange(direction=RIGHT))
    
        
#python -m manim -psqh illustration.py Networks