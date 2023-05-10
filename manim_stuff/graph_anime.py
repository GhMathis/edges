# -*- coding: utf-8 -*-
"""
Created on Tue May  9 17:52:49 2023

@author: Mathis
"""

from manim import *
import numpy as np

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
        
    def construct(self):
        adj_arr = np.array([[0, 1, 1, 1, 1],
                            [1, 0, 1, 0, 1],
                            [0, 1, 0, 1, 0],
                            [1, 0, 1, 0, 1],
                            [1, 1, 0, 1, 0]])
       
        ##### Setup the data for visualisations
        self.adj_to_list(adj_arr)
        
        
        meta_ntw = Graph(vertices = self.vertices,
                         edges =self.edges )
        
        self.play(LaggedStart(Create(meta_ntw),  run_time=3))
        self.wait(2)
        
        
#python -m manim -pql graph_anime.py Networks