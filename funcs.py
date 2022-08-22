##########
#
# funcs.py
#
#
# Author: Clinton H. Durney
# Email: clinton.durney@jic.ac.uk 
#
# Last Edit: 08/22/22
##########

import networkx as nx
from scipy.spatial import distance
import numpy as np
import globals as const
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pdb
import sys

###############
def tissue_2d():

    def gen_nodes(ori):
        nodes = [[ori[0] + r*np.cos(n*np.pi/3), ori[1] + r*np.sin(n*np.pi/3)] for n in range(0,6)]
        return nodes

    def add_nodes(nodes, i):
        pos = nx.get_node_attributes(G,'pos')
        cen_index = i-1
        if i < 1000:
            centers.append(cen_index)
        AS_boundary = []
        spokes = []
        for node in nodes:
            add_node = True
            for existing_node in pos:
                if distance.euclidean(pos[existing_node],node) < 10**(-7):
                    add_node = False
                    AS_boundary.append(existing_node)
                    spokes.append((cen_index,existing_node))
                    break

            if add_node == True:
                G.add_node(i,pos=node)
                i += 1
                AS_boundary.append(i-1)
                spokes.append((cen_index,i-1))

        return AS_boundary, spokes, i

    def add_spokes_edges(spokes, boundary):
        boundary.append(boundary[0])
        G.add_edges_from(spokes, l_rest=const.l_apical, myosin=0)    
        G.add_path(boundary, l_rest = const.l_apical, myosin=0, color='#808080')

        return

    G = nx.Graph()

    r = const.l_apical              # initial spoke length
    num_cells = 2*const.hex-1          # number of cells in center row

    centers = []
    
    # Apical Nodes
    i = 0
    # Center cell set up
    origin = [0.0,0.0]
    G.add_node(i,pos=origin)
    i += 1

    nodes = gen_nodes(origin)
    AS_boundary, spokes, i = add_nodes(nodes,i)
    add_spokes_edges(spokes, AS_boundary)

    for index in range(1,int((num_cells - 1)/2.)+1):
        # # Step Up
        origin = [0, np.sqrt(3)*r*index]
        G.add_node(i,pos=origin)
        i += 1

        nodes = gen_nodes(origin)
        AS_boundary, spokes, i = add_nodes(nodes,i)
        add_spokes_edges(spokes, AS_boundary)

        # # # Step down
        origin = [0, -np.sqrt(3)*r*index]
        G.add_node(i,pos=origin)
        i += 1

        nodes = gen_nodes(origin)
        AS_boundary, spokes, i = add_nodes(nodes,i)
        add_spokes_edges(spokes, AS_boundary)

    for index in range(1,const.hex):  
        if (num_cells - index) % 2 == 0:
            for j in range(1,(num_cells-index),2):
                origin = [(3/2.)*r*index,(np.sqrt(3)/2.)*r*j]
                G.add_node(i,pos=origin)
                i += 1

                nodes = gen_nodes(origin)
                AS_boundary, spokes, i = add_nodes(nodes,i)
                add_spokes_edges(spokes, AS_boundary)

                origin = [(3/2.)*r*index,(-np.sqrt(3)/2.)*r*j]
                G.add_node(i,pos=origin)
                i += 1

                nodes = gen_nodes(origin)
                AS_boundary, spokes, i = add_nodes(nodes,i)
                add_spokes_edges(spokes, AS_boundary)

            # Step Left

                origin = [-(3/2.)*r*index,(np.sqrt(3)/2.)*r*j]
                G.add_node(i,pos=origin)
                i += 1

                nodes = gen_nodes(origin)
                AS_boundary, spokes, i = add_nodes(nodes,i)
                add_spokes_edges(spokes, AS_boundary)

                origin = [-(3/2.)*r*index,(-np.sqrt(3)/2.)*r*j]
                G.add_node(i,pos=origin)
                i += 1

                nodes = gen_nodes(origin)
                AS_boundary, spokes, i = add_nodes(nodes,i)
                add_spokes_edges(spokes, AS_boundary)

        else:
            for j in range(0,(num_cells-index),2):
                origin = [3*(1/2.)*r*index, (np.sqrt(3)/2.)*r*j]
                G.add_node(i,pos=origin)
                i += 1

                nodes = gen_nodes(origin)
                AS_boundary, spokes, i = add_nodes(nodes,i)
                add_spokes_edges(spokes, AS_boundary)
                
                if j != 0:
                    origin = [3*(1/2.)*r*index, -(np.sqrt(3)/2.)*r*j]
                    G.add_node(i,pos=origin)
                    i += 1

                    nodes = gen_nodes(origin)
                    AS_boundary, spokes, i = add_nodes(nodes,i)
                    add_spokes_edges(spokes, AS_boundary)

                # Step Left
                origin = [-3*(1/2.)*r*index, (np.sqrt(3)/2.)*r*j]
                G.add_node(i,pos=origin)
                i += 1

                nodes = gen_nodes(origin)
                AS_boundary, spokes, i = add_nodes(nodes,i)
                add_spokes_edges(spokes, AS_boundary)
                
                if j != 0:
                    origin = [-3*(1/2.)*r*index, -(np.sqrt(3)/2.)*r*j]
                    G.add_node(i,pos=origin)
                    i += 1

                    nodes = gen_nodes(origin)
                    AS_boundary, spokes, i = add_nodes(nodes,i)
                    add_spokes_edges(spokes, AS_boundary)
   
    circum_sorted = []
    pos = nx.get_node_attributes(G,'pos') 
    xy = [[pos[n][0],pos[n][1]] for n in range(0,i)]
    for center in centers:
        a, b = sort_corners(list(G.neighbors(center)),xy[center],xy)
        circum_sorted.append(np.asarray([b[n][0] for n in range(len(b))]))
    circum_sorted = np.array(circum_sorted)

    belt = []
    for node in G.nodes():
        if len(list(G.neighbors(node))) < 6:
            belt.append(node)
    xy_belt = [xy[n] for n in belt]
    a,b = sort_corners(belt, xy[0], xy)
    belt = np.array([b[n][0] for n in range(len(b))])
    
    # make swaps for the "corner" nodes that the angle doesn't account for
    for n in range(1,len(belt)-1):
        if G.has_edge(belt[n-1],belt[n]) == False:
            belt[n], belt[n+1] = belt[n+1], belt[n]

    triangles = []
    for node in G.nodes():
        if node not in belt:
            if node in centers:
                out1, out2 = sort_corners(list(G.neighbors(node)),pos[node],pos)
                neighbors = [out2[k][0] for k in range(0,len(out2))]
                alpha_beta = [[[node,neighbors[k-1],neighbors[k-2]],[node, neighbors[k],neighbors[k-1]]] for k in range(0,6)]

                for entry in alpha_beta:
                    triangles.append(entry)
            else: # node not a center, so that I don't double count pairs, only keep those that cross a cell edge
                out1, out2 = sort_corners(list(G.neighbors(node)),pos[node],pos)
                neighbors = [out2[k][0] for k in range(0,len(out2))]
	            
                for k in range(0,6):
                    alpha = [node,neighbors[k-1],neighbors[k-2]]
                    beta = [node,neighbors[k],neighbors[k-1]]
                    
                    if set(alpha) & set(centers) != set(beta) & set(centers):
                        triangles.append([alpha,beta])

    print("Apical nodes added correctly.")
    print("Number of apical nodes are", i)
    
    return G, centers, i, circum_sorted, belt, triangles
###############

def vector(A,B):
        
    return [(B[0]-A[0]), (B[1]-A[1])] 

def unit_vector(A,B):
    # Calculate the unit vector from A to B in 3D

    dist = distance.euclidean(A,B)

    if dist < 10e-15:
        dist = 1.0

    return [(B[0]-A[0])/dist,(B[1]-A[1])/dist]
###############

def d_pos(position,force,dt):
    # Calculate the new vertex position using forward euler
    # input: Current position, force, and dt the time step
    # output: Updated position of the node.  

    x_new = position[0] + (dt/const.eta)*force[0]
    y_new = position[1] + (dt/const.eta)*force[1]
    
    return [x_new,y_new]
###############

def elastic_force(l,l0,muu):
    # Calculate the magnitude of the force obeying Hooke's Law

    frce = muu*(l-l0) 

    return frce 
###############

def signed_angle(v1,v2):
    theta = np.arctan2(v2[1],v2[0]) - np.arctan2(v1[1],v1[0])
    if theta > np.pi:
        theta -= 2*np.pi
    elif theta <= -np.pi:
        theta += 2*np.pi
    return theta
###############

def get_points(G, q, pos):
    # get node numbers associated with a given center
    # inputs:   G: networkx graph
    #           q: number of center node (apical only)
    #           pos: position of nodes 
    # returns:  pts: list of positions that are associated with that center 

    pts = [q] + list(G.neighbors(q))
    pts = [pos[n] for n in pts]

    return pts 

def sort_corners(corners,center_pos,pos_nodes):

    corn_sort = [(corners[0],0)]
    u = unit_vector(center_pos,pos_nodes[corners[0]])
    for i in range(1,len(corners)):
        v = unit_vector(center_pos,pos_nodes[corners[i]])
        dot = np.dot(u,v)
        det = np.linalg.det([u,v])
        angle = np.arctan2(det,dot)
        corn_sort.append((corners[i],angle))
        corn_sort = sorted(corn_sort, key=lambda tup: tup[1])
        corn2 = [pos_nodes[entry[0]] for entry in corn_sort]
    
    return corn2, corn_sort

def CW(x):

    return [x[1],-x[0]]

def CCW(x):
    
    return [-x[1],x[0]]

def innie(theta, u, v):

    if theta >= const.alpha:
        udir = CCW(u)
        vdir = CW(v)
    else:
        udir = CW(u)
        vdir = CCW(v)
    
    return udir, vdir

def outie(theta, u, v):

    udir = CCW(u)
    vdir = CW(v)

    return udir, vdir




