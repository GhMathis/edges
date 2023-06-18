# -*- coding: utf-8 -*-
"""
Created on Tue May  9 17:52:49 2023

@author: Mathis
"""

from manim import *
from manim.utils.color import Colors
import numpy as np
import networkx as nx
import math as m
from scipy.linalg import expm
class Networks(Scene):

    config.background_color = WHITE
    def adj_to_list(self,adj_arr):
        """Converte an adjacency matrix in the form of list of vertices and edges
        vertices = list(v1, v2, v3,...)
        edges = list((v1, v2), (v3, v1),...) """
        
        ##### vertices conversion
        vertices = np.arange(adj_arr.shape[0]).tolist()
        ID = np.triu_indices(adj_arr.shape[0])
        
        ##### Edges conversion
        adj_arr[ID] = 0
        edges_temp = np.where(adj_arr == 1)
        edges = [(e1,e2) for e1,e2 in zip(edges_temp[0],edges_temp[1])]
        return([vertices, edges])
    
    def egde_color(self,decompo_arr):
        """Color depending on the communicability group"""
        ID = np.triu_indices(decompo_arr.shape[0],1)
        decompo_arr[ID] = 0
        pos_edges =  np.where(decompo_arr > 0)
        neg_edges =  np.where(decompo_arr < 0)
        egde_col_pos = {(e1,e2):{"stroke_color": Colors.blue.value} for e1,e2 in zip(pos_edges[0], pos_edges[1])}
        egde_col_pos2 = {(e1,e2):{"stroke_color": Colors.red.value} for e1,e2 in zip(neg_edges[0], neg_edges[1])}
        egde_col_pos = egde_col_pos| egde_col_pos2
        
        return(egde_col_pos)
    def egde_black(self,edges, color = BLACK):
        """Color all in one (BLACK default)"""
        if type(color) == list:
            egde_col = {edges_tuple:{"stroke_color": c}
                            for edges_tuple,c in zip(edgesc,color)}
        else:
            egde_col = {edges_tuple:{"stroke_color": color}
                            for edges_tuple in edges}
        return(egde_col)
    
    def vertex_color(self,vertex,fill=["#A4E8FF"] ):
        """Color all in one (BLACK default)"""
        if type(fill) == list :
            
            verte_col = {v:{"fill_color": f} for v, f in zip(vertex,fill)}
        else:   
            verte_col = {v:{"fill_color": fill} for v in vertex}
        
        return(verte_col)
    def verte_color(self,eign_vec):
        """Color depending on the communicability group"""
        pos_verte =  np.where(eign_vec > 0)
        neg_verte =  np.where(eign_vec < 0)
        
        verte_col_pos = {v:{"fill_color": Colors.blue.value} for v in pos_verte[0]}
        verte_col_pos2 = {v:{"fill_color": Colors.red.value} for v in neg_verte[0]}
        verte_col_pos = verte_col_pos| verte_col_pos2
        
        return(verte_col_pos)
    


        
    def construct(self):

        def all_arrows(lright, lleft, k, color):
            """Arrow Vgroup for anniimation"""
            arrows = []
            
            for r,l,c in zip(lright,lleft,color):
                arrows.append(Line(start = r + tracker.get_value()*(l-r)*
                               (np.sqrt(tracker.get_value())),
                               end = r + tracker.get_value()*(l-r)*
                               (2-tracker.get_value()),
                               color = c, buff = 0, tip_length = 0.3 *k).add_tip()
                              )
                
            return VGroup(*arrows)
        
        
        adj_arr = np.array([[0, 0, 0, 1, 1, 0, 0],
                            [0, 0, 0, 0, 1, 1, 1],
                            [0, 0, 0, 0, 0, 1, 1],
                            [1, 0, 0, 0, 0, 0, 0],
                            [1, 1, 0, 0, 0, 0, 0],
                            [0, 1, 1, 0, 0, 0, 0],
                            [0, 1, 1, 0, 0, 0, 0],])
        
        spectra = np.linalg.eig(adj_arr)
        ##### Setup the data for visualisations
        vertices1, edges1 = self.adj_to_list(adj_arr)
        
        
        G_decomposition = lambda val, vec: np.dot(vec.T, vec)*np.exp(val)
        list_g_decomp = [G_decomposition(val, vec[np.newaxis]) for val, vec
                 in zip(spectra[0], spectra[1])]
        list_g_decomp[0] == list_g_decomp[0].T
        
        ##### matrix 1
        G_mat = []
        for i in range(5):
            if i != 0:
                G_mat.append(np.dot(np.dot(spectra[1], np.diag(spectra[0]**i/m.factorial(i))),
                       np.transpose(spectra[1]))+G_mat[-1])
            else:
                G_mat.append(np.dot(np.dot(spectra[1], np.diag(spectra[0]**i/m.factorial(i))),
                        np.transpose(spectra[1])))
        object_mat1 = []
        for i in range(2):
            G_mat[i] = np.round(G_mat[i],2)
            object_mat1.append(Matrix(list(G_mat[i]),).scale(0.5))
        print(object_mat1[1])
        background= Rectangle().scale(1.6)
        background.set_fill(opacity=.5)
        background.set_color(BLACK)
            
        ##### matrix 2
        adj_arr2 = np.array([[0, 0, 0, 1, 1, 0, 0],
                            [0, 0, 0, 0, 1, 0, 1],
                            [0, 0, 0, 0, 0, 1, 1],
                            [1, 0, 0, 0, 0, 0, 0],
                            [1, 1, 0, 0, 0, 0, 0],
                            [0, 0, 1, 0, 0, 0, 0],
                            [0, 1, 1, 0, 0, 0, 0],])
        
        spectra2 = np.linalg.eig(adj_arr2)
        vertices2, edges2 = self.adj_to_list(adj_arr2)
        
        # G_mat2 = []
        # for i in range(5):
        #     if i != 0:
        #         G_mat2.append(np.dot(np.dot(spectra2[1], np.diag(spectra2[0]**i/m.factorial(i))),
        #                np.transpose(spectra2[1]))+G_mat2[-1])
        #     else:
        #         G_mat2.append(np.dot(np.dot(spectra2[1], np.diag(spectra2[0]**i/m.factorial(i))),
        #                 np.transpose(spectra2[1])))
        # for i in range(5):
        #     G_mat2[i] = np.abs(np.round(G_mat2[i],2))
            
        ##### Object creation
        
        
        ### network
        lyt = nx.bipartite_layout(vertices1[:3], vertices1[3:])
        G1 = nx.Graph()
        G1.add_nodes_from(vertices1)
        G1.add_edges_from(edges1)
        meta_ntw1 = Graph(list(G1.nodes), list(G1.edges),
                         labels =True,
                         layout = "partite",
                         edge_config = self.egde_black(list(G1.edges)),
                         vertex_config = self.vertex_color(list(G1.nodes),
                         fill=["#808080", "#808080", "808080",
                         "#FB8072", "808080", "#FDB462", "808080"]),
                         partitions = [[3,4,5,6]]) #for g_decomp
                    #in list_g_decomp[:5]]
        meta_ntw1bis = Graph(list(G1.nodes), list(G1.edges),
                         labels =True,
                         layout = "partite",
                         edge_config = self.egde_black(list(G1.edges)),
                         vertex_config = self.vertex_color(list(G1.nodes),
                         fill=["#808080", "#808080", "#808080",
                         "#808080", "#808080", "#FDB462", "#B3DE69"]),
                         partitions = [[3,4,5,6]]) #for g_decomp
        meta_ntw1bis2 = Graph(list(G1.nodes), list(G1.edges),
                         labels =True,
                         layout = "partite",
                         edge_config = self.egde_black(list(G1.edges)),
                         vertex_config = self.vertex_color(list(G1.nodes),
                         fill=["#808080", "#808080", "#BEBADA",
                         "#808080", "#808080", "#808080", "#B3DE69"]),
                         partitions = [[3,4,5,6]])
        #["#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5"]  
        lyt = nx.bipartite_layout(vertices2[:3], vertices2[3:])
        G2 = nx.Graph()
        G2.add_nodes_from(vertices2)
        G2.add_edges_from(edges2)
        meta_ntw2= Graph(list(G2.nodes), list(G2.edges),
                         labels =True,
                         layout = "partite",
                         edge_config = self.egde_black(list(G2.edges)),
                         vertex_config = self.vertex_color(list(G2.nodes),
                         fill=["#808080", "#808080", "#808080",
                         "#808080", "#808080", "#FDB462", "#B3DE69"]),
                         partitions = [[3,4,5,6]])#for g_decomp 
        graphs  = [meta_ntw1,meta_ntw2]
        r1 = VGroup(*graphs).arrange_in_grid(1, buff=1)
        #r2 = VGroup(meta_ntw2).arrange(direction=DOWN)

        #path_list = VGroup(*path_list)
        edge_node_dict = {}

        #edge_node_dict = {f"{n}":[r1[e+7] for e in meta_ntw1[n].keys()] for n in  range(7)}
       
        for vertex in meta_ntw1.vertices:
            connected_edges = []
            for edge in meta_ntw1.edges:
                if vertex in edge:
                    connected_edges.append(edge)
      
        #move_sequence = AnimationGroup(*move_animations)
        
        ####### Animation #######
        ### setup
     
       
        
        ### Start
        # self.play(LaggedStart(VGroup(r1).shift(4*LEFT),
        #                        run_time=2))
        max_step = 2
        # ##### 1
        # self.add(meta_ntw1.shift(4*LEFT))
        # self.wait(1)
        
        # start_temp = [0,1,2,3,4,5,6]
        # move_animations= {}
        # col_vertex = [v.fill_color for
        #               v in list(meta_ntw1.vertices.values())]
        
        # start_color_temp = []
        # for i in start_temp:
        #     start_color_temp.append(col_vertex[i])
            
        # # start_color_temp = [(col_vertex[e[0]], col_vertex[e[1]]) for
        # #                e in list(meta_ntw1.edges.keys())]
        
        # for step in range(max_step):
        #     tracker = ValueTracker(0.009)
            
        #     start_color = start_color_temp
        #     start_color_temp = []
        #     arrow_color= []
            
        #     start = start_temp
        #     #print(len(start), len(start_color))
        #     start_temp = []
            
        #     lleft = []
        #     lright = []
            
        #     lines = list(meta_ntw1.edges.values())
        #     edge_id = list(meta_ntw1.edges.keys())

        #     for i, c in zip(start, start_color):
        #         connect_edges = []
        #         for e, l in zip(edge_id, lines):
        #             if i == e[0]:
        #                 #print(["v:", v, " start : " ,start, "step : ", step ])    
        #                 lleft.append(l.get_end())
        #                 lright.append(l.get_start())
        #                 start_color_temp.append(c)
        #                 #arrow_color.append(c[0])
        #                 start_temp.append(e[1])
        #             if  i == e[1]:
        #                     #print(["v:", v, " start : " ,start, "step : ", step ])
        #                 lleft.append(l.get_start())
        #                 lright.append(l.get_end())
        #                 start_color_temp.append(c)
        #                 #arrow_color.append(c[1])
        #                 start_temp.append(e[0])
           
        #     draw_arrows = always_redraw(
        #     lambda: all_arrows(lright, lleft, 1/m.factorial(step+1),
        #                        color = start_color_temp))
            
        #     #print(start_color_temp)
        #     #print(start_temp)
        #     self.play(Create(draw_arrows))
        #     self.play(tracker.animate.set_value(1), run_time=2, rate_func = linear)
            
        # self.wait(1)
            
        
        
        
        # ##### 2
        self.play(Create(meta_ntw1bis.shift(LEFT)))#.next_to(meta_ntw1, RIGHT, buff = 1)))
        
        start_temp = [5]
        move_animations= {}
        col_vertex = [v.fill_color for
                      v in list(meta_ntw1bis.vertices.values())]
        
        start_color_temp = []
        for i in start_temp:
            start_color_temp.append(col_vertex[i])
            
        # start_color_temp = [(col_vertex[e[0]], col_vertex[e[1]]) for
        #                e in list(meta_ntw1bis.edges.keys())]
        
        for step in range(max_step):
            tracker = ValueTracker(0.009)
            
            start_color = start_color_temp
            start_color_temp = []
            arrow_color= []
            
            start = start_temp
            #print(len(start), len(start_color))
            start_temp = []
            
            lleft = []
            lright = []
            
            lines = list(meta_ntw1bis.edges.values())
            edge_id = list(meta_ntw1bis.edges.keys())

            for i, c in zip(start, start_color):
                connect_edges = []
                for e, l in zip(edge_id, lines):
                    if i == e[0]:
                        #print(["v:", v, " start : " ,start, "step : ", step ])    
                        lleft.append(l.get_end())
                        lright.append(l.get_start())
                        start_color_temp.append(c)
                        #arrow_color.append(c[0])
                        start_temp.append(e[1])
                    if  i == e[1]:
                            #print(["v:", v, " start : " ,start, "step : ", step ])
                        lleft.append(l.get_start())
                        lright.append(l.get_end())
                        start_color_temp.append(c)
                        #arrow_color.append(c[1])
                        start_temp.append(e[0])
           
            draw_arrows = always_redraw(
            lambda: all_arrows(lright, lleft, 1/m.factorial(step+1),
                                color = start_color_temp))
            
            #print(start_color_temp)
            #print(start_temp)
            self.play(Create(draw_arrows))
            self.play(tracker.animate.set_value(1), run_time=4, rate_func = smooth)
            
        self.wait(1)
        # ##### 3
        self.play(Create(meta_ntw2.next_to(meta_ntw1bis, RIGHT, buff = 1)))
            
        start_temp = [5]
        move_animations= {}
        col_vertex = [v.fill_color for
                      v in list(meta_ntw2.vertices.values())]
        
        start_color_temp = []
        for i in start_temp:
            start_color_temp.append(col_vertex[i])
            
        # start_color_temp = [(col_vertex[e[0]], col_vertex[e[1]]) for
        #                e in list(meta_ntw2.edges.keys())]
        
        for step in range(max_step):
            tracker = ValueTracker(0.009)
            
            start_color = start_color_temp
            start_color_temp = []
            arrow_color= []
            
            start = start_temp
            #print(len(start), len(start_color))
            start_temp = []
            
            lleft = []
            lright = []
            
            lines = list(meta_ntw2.edges.values())
            edge_id = list(meta_ntw2.edges.keys())

            for i, c in zip(start, start_color):
                connect_edges = []
                for e, l in zip(edge_id, lines):
                    if i == e[0]:
                        #print(["v:", v, " start : " ,start, "step : ", step ])    
                        lleft.append(l.get_end())
                        lright.append(l.get_start())
                        start_color_temp.append(c)
                        #arrow_color.append(c[0])
                        start_temp.append(e[1])
                    if  i == e[1]:
                            #print(["v:", v, " start : " ,start, "step : ", step ])
                        lleft.append(l.get_start())
                        lright.append(l.get_end())
                        start_color_temp.append(c)
                        #arrow_color.append(c[1])
                        start_temp.append(e[0])
            
            draw_arrows = always_redraw(
            lambda: all_arrows(lright, lleft, 1/m.factorial(step+1),
                                color = start_color_temp))
            
            #print(start_color_temp)
            #print(start_temp)
            self.play(Create(draw_arrows))
            self.play(tracker.animate.set_value(1), run_time=5, rate_func = smooth)
            
        self.wait(1)
        
        
        ##### 4

        # self.play(Create(meta_ntw1bis2))#.next_to(meta_ntw2, RIGHT, buff = 1)))
        
        # start_temp = [2,6]
        # move_animations= {}
        # col_vertex = [v.fill_color for
        #               v in list(meta_ntw1bis2.vertices.values())]
        
        # start_color_temp = []
        # for i in start_temp:
        #     start_color_temp.append(col_vertex[i])
            
        # # start_color_temp = [(col_vertex[e[0]], col_vertex[e[1]]) for
        # #                e in list(meta_ntw1bis2.edges.keys())]
        # arrow_list = []
        # for step in range(max_step):
        #     tracker = ValueTracker(0.009)
            
        #     start_color = start_color_temp
        #     start_color_temp = []
        #     arrow_color= []
            
        #     start = start_temp
        #     #print(len(start), len(start_color))
        #     start_temp = []
            
        #     lleft = []
        #     lright = []
            
        #     lines = list(meta_ntw1bis2.edges.values())
        #     edge_id = list(meta_ntw1bis2.edges.keys())
        #     if step == 0:
        #         lleft.append(lines[6].get_end())
        #         lright.append(lines[6].get_start())
        #         start_color_temp.append(start_color[0])
        #         start_temp.append(edge_id[6][1])
                
        #         lleft.append(lines[6].get_start())
        #         lright.append(lines[6].get_end())
        #         start_color_temp.append(start_color[1])
        #         #arrow_color.append(c[1])
        #         start_temp.append(edge_id[6][0])


        #     else:

        #         for i, c in zip(start, start_color):
        #             connect_edges = []
        #             for e, l in zip(edge_id, lines):
        #                 if i == e[0]:
        #                     #print(["v:", v, " start : " ,start, "step : ", step ])    
        #                     lleft.append(l.get_end())
        #                     lright.append(l.get_start())
        #                     start_color_temp.append(c)
        #                     #arrow_color.append(c[0])
        #                     start_temp.append(e[1])
        #                 if  i == e[1]:
        #                         #print(["v:", v, " start : " ,start, "step : ", step ])
        #                     lleft.append(l.get_start())
        #                     lright.append(l.get_end())
        #                     start_color_temp.append(c)
        #                     #arrow_color.append(c[1])
        #                     start_temp.append(e[0])
        #     arrow_list.append([lright,lleft, start_color_temp])
        
        
        # for step in range(max_step):   
        #     tracker = ValueTracker(0.009)
        #     draw_arrows = always_redraw(
        #     lambda: all_arrows(arrow_list[step][0], arrow_list[step][1], 1/(step+1),
        #                         color = arrow_list[step][2]))
        #     print(arrow_list[step])   
        #     #print(start_color_temp)
        #     #print(start_temp)
        #     self.play(Create(draw_arrows))
        #     self.play(tracker.animate.set_value(1), run_time=2, rate_func = linear)
            
            # self.wait(1)
            # move_sequence = AnimationGroup(*move_animations)
            # move_animations[f"{step}"]= VGroup(*move_animation)
            # self.play(move_sequence, run_time=10, rate_func=linear)
        
# python -m manim -pql graph_anime.py Networks     
        
# [[1,6,2,5],[1,4,1,5],[1,6,1,5], [1,5,1,5], [1,5,2,5],
#  [5,2,6,1],[5,1,4,1], [5,1,6,1], [5,1,5,1], [5,2,4]]