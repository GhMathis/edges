# -*- coding: utf-8 -*-
"""
Created on Tue May  9 17:52:49 2023

@author: Mathis
"""

from manim import *
from manim.utils.color import Colors
import numpy as np
import networkx as nx
import math

class Networks(Scene):
    
    def adj_to_list(self,adj_arr):
        """Converte an adjacency matrix in the form of list of vertices and edges
        vertices = list(v1, v2, v3,...)
        edges = list((v1, v2), (v3, v1),...) """
        
        ##### vertices conversion
        self.vertices = np.arange(adj_arr.shape[0]).tolist()
        ID = np.triu_indices(adj_arr.shape[0])
        
        ##### Edges conversion
        adj_arr[ID] = 0
        edges_temp = np.where(adj_arr == 1)
        self.edges = [(e1,e2) for e1,e2 in zip(edges_temp[0],edges_temp[1])]
    
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

        def arrowupdate_toright(arrow, l, t):
            start = l.get_end()
            end = l.get_start()

            p = start+t*(end-start)*(2-t)
            p2 =  start+t*(end-start)*(np.sqrt(t))
            #start = l.get_end()

            arrow.put_start_and_end_on(p2, p)
        def trackfunc_toright(arrow,l):
            t = tracker.get_value()
            arrowupdate_toright(arrow, l, t)

        def arrowupdate_toleft(arrow, l, t):
             start = l.get_start()
             end = l.get_end()

             p = start+t*(end-start)*(2-t)
             p2 =  start+t*(end-start)*(np.sqrt(t))
             #start = l.get_end()

             arrow.put_start_and_end_on(p2, p)
             

        def trackfunc_toleft(arrow,l):
            t = tracker.get_value()
            arrowupdate_toleft(arrow, l, t)    
        
        adj_arr = np.array([[0, 0, 0, 1, 1, 0, 0],
                            [0, 0, 0, 0, 1, 1, 1],
                            [0, 0, 0, 0, 0, 1, 1],
                            [1, 0, 0, 0, 0, 0, 0],
                            [1, 1, 0, 0, 0, 0, 0],
                            [0, 1, 1, 0, 0, 0, 0],
                            [0, 1, 1, 0, 0, 0, 0],])
        
        spectra = np.linalg.eig(adj_arr)
        ##### Setup the data for visualisations
        self.adj_to_list(adj_arr)
        
        
        G_decomposition = lambda val, vec: np.dot(vec.T, vec)*np.exp(val)
        list_g_decomp = [G_decomposition(val, vec[np.newaxis]) for val, vec
                 in zip(spectra[0], spectra[1])]
        list_g_decomp[0] == list_g_decomp[0].T

        ##### Object creation
        lyt = nx.bipartite_layout(self.vertices[:3], self.vertices[3:])
        G = nx.Graph()
        G.add_nodes_from(self.vertices)
        G.add_edges_from(self.edges)
        meta_ntw = Graph(list(G.nodes), list(G.edges),
                         labels =True,
                         layout = "partite",
                         #edge_config = self.egde_color(g_decomp),
                         partitions = [[3,4,5,6]]) #for g_decomp
                    #in list_g_decomp[:5]]
        
        
        r1 = VGroup(*meta_ntw)
        ###
        d1 = [Dot() for _ in range(len(r1[7:])*2)]
        d1 = VGroup(*d1)
        #path_list = VGroup(*path_list)
        edge_node_dict = {}

        #edge_node_dict = {f"{n}":[r1[e+7] for e in meta_ntw[n].keys()] for n in  range(7)}
       
        for vertex in meta_ntw.vertices:
            connected_edges = []
            for edge in meta_ntw.edges:
                if vertex in edge:
                    connected_edges.append(edge)
      
        #move_sequence = AnimationGroup(*move_animations)
        start = [1]
        start_temp = [1]
        tracker = ValueTracker(0.009)
        # print(meta_ntw.edges)
        # l = meta_ntw.edges[(0,3)]
        
        # lend= l.get_end()
        # lstart = l.get_start()
        # p = lstart+0.01*(lend-lstart)
        # A = Arrow(start = lstart, end = p, color = GREEN, buff = 0,
        #           max_stroke_width_to_length_ratio = 20,
        #           max_tip_length_to_length_ratio=20)
        # move_animation = A.add_updater(lambda arrow : trackfunc(arrow,l)).update()

        move_animations= {}
        print(meta_ntw.edges.items())
        def all_arrows(lright, lend, k):
            arrows = []
            for r,l in zip(lright,lleft):
               arrows.append(Line(start = r + tracker.get_value()*(l-r)*
                              (np.sqrt(tracker.get_value())),
                              end = r + tracker.get_value()*(l-r)*
                              (2-tracker.get_value()),
                              color = GREEN, buff = 0, tip_length = 0.3 *k
                          # max_stroke_width_to_length_ratio = 20,
                          # max_tip_length_to_length_ratio=20*k,
                          ).add_tip()
                             )
            return VGroup(*arrows)
        
        self.play(LaggedStart(Create(r1),
                              run_time=3))
        self.wait(1)
    
        for step in range(4):
            tracker = ValueTracker(0.009)
            start = start_temp
            start_temp = []
            move_animation = []
            lleft = []
            lright = []
            for v,l in meta_ntw.edges.items():
            
                if v[0] in start and v[0] in [0,1,2]:

                    lleft.append(l.get_end())
                    lright.append(l.get_start())
                
                elif v[1] in start and v[1] in [3,4,5,6]:
                    print(["v:", v, " start : " ,start, "step : ", step ])
                    lleft.append(l.get_start())
                    lright.append(l.get_end())

                draw_arrows = always_redraw(
                    lambda: all_arrows(lright, lleft, 1/math.factorial(step+1)))

                print([" start : " ,start, "step : ", step ])
                if v[0] in start or v[1] in start:
                    if v[0] not in start:
                        start_temp.append(v[0])
                    if v[1] not in start:
                        start_temp.append(v[1])
            self.play(Create(draw_arrows))
            self.play(tracker.animate.set_value(1), run_time=2, rate_func = linear)
            #move_sequence = AnimationGroup(*move_animations)
            #move_animations[f"{step}"]= VGroup(*move_animation)
            #self.play(move_sequence, run_time=10, rate_func=linear)
        
# python -m manim -pql graph_anime.py Networks     
        
