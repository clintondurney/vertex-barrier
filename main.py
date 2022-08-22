import networkx as nx
import numpy as np
from scipy.spatial import distance
import csv
import globals as const
from funcs import *

# Constants for simulation
dt = const.dt
#var_dt = True

# dimensions of the cell 
l_apical = const.l_apical 

# Set the arcs list
inner_arc = const.inner_arc
outer_arc = const.outer_arc

# mechanical parameters
l0_apical = l_apical
mu_apical = const.mu_apical         
myo_beta = const.myo_beta 
eta = const.eta 

# initialize the tissue
G, centers, num_api_nodes, circum_sorted, belt, triangles = tissue_2d()
pit_centers = const.pit_centers 

# Starting from t=0
t = 0 

# Starting from t != 0.  Read previous .pickle file.
# Useful if doing a "manual" shortening of the rest length to find a ground state configuration.
#t = 1000
#G = nx.read_gpickle('/home/cdurney/dina-models/2d-vertex/outer-stiff/arc-short/t1000.pickle')

# t=initial nx Graph in pickled form for plotting later
print(t) 
file_name = 't' + str(int(t)) 
nx.write_gpickle(G,file_name + '.pickle')
np.save(file_name,circum_sorted) 

while t <= const.t_final:
    
    # increment t by dt
    # initialize force_dict back to zeros
    t = round(t+dt,1)
    print(dt, t) 
    pos = nx.get_node_attributes(G,'pos')
    force_dict = {new_list: np.zeros(2,dtype=float) for new_list in G.nodes()} 

    # Update myosin on "pit cells" i.e., contraction in the center of the tissue caused by compression
    if t == const.t_pit: 
        for node in pit_centers: 
            if node == 0:
                myo = 1.39*const.pit_strength
            for neighbor in G.neighbors(node): 
                G[node][neighbor]['myosin'] = const.pit_strength 
        print("Pit is established")

    # update myosin on inner arc 
    if t == const.t_1:
        for i in range(0,len(inner_arc)):
            G[inner_arc[i-1]][inner_arc[i]]['myosin'] = const.belt_strength     
        print("Inner arc established")

    # update myosin on outer arc -- not needed for this work, but kept for adapting this code
#    if t == const.t_2:
#        for i in range(0,len(outer_arc)):
#            G[outer_arc[i-1]][outer_arc[i]]['myosin'] = const.arc_strength     
#        print("Outer arc established")

#    # update myosin on belt -- not needed for this work, but kept for adapting this code
#    if t == const.t_belt:
#        for i in range(0,len(belt)):
#            G[belt[i-1]][belt[i]]['myosin'] = const.belt_strength     
#        print("Belt established") 

    for node in G.nodes(): 
        # update force on each node  
        force = [0.0,0.0]
    
        # Elastic forces due to the cytoskeleton 
        for neighbor in G.neighbors(node):
            a = pos[node]
            b = pos[neighbor]
            
            dist = distance.euclidean(a,b)
            direction = unit_vector(a,b)
            
            if (node in outer_arc) and (neighbor in outer_arc):
                magnitude = elastic_force(dist, const.ring_length, const.mu_stiff) 
            else:
                magnitude = elastic_force(dist, G[node][neighbor]['l_rest'], mu_apical) 
            force = np.sum([force,magnitude*np.array(direction)],axis=0)
            
            # Force due to myosin
            magnitude = myo_beta*G[node][neighbor]['myosin']
            force = np.sum([force, magnitude*np.array(direction)],axis=0)

        force_dict[node] = np.add(force_dict[node], force) 

    # Implement bending energy on outer arc -- i.e., the "stiff boundary"
    for i in range(-1,len(outer_arc)-1):
        a = pos[outer_arc[i-1]]
        b = pos[outer_arc[i]]
        c = pos[outer_arc[i+1]]
        
        u = unit_vector(b,c)
        v = unit_vector(b,a)
        
        theta = signed_angle(u,v)
        if theta >= 0:
            udir, vdir = innie(theta, u, v)
        else:
            udir, vdir = outie(theta, u, v) 

        magnitude = (1/4)*const.k*np.abs(np.abs(theta)-const.alpha)         
        # factor of 1/4, because they need to be split across the nodes. 
        # (1/4)+(1/4) = (1/2) of the bending force to each branch of the junction 
         
        # i-1 node goes in v-dir 
        force = np.sum([[0.0,0.0], magnitude*np.array(vdir)],axis=0)
        force_dict[outer_arc[i-1]] = np.add(force_dict[outer_arc[i-1]], force) 
        
        # i node goes in sum of u+v/|u+v|
        sum_uv = [u[0]+v[0],u[1]+v[1]]
        sum_uv = sum_uv/np.linalg.norm(sum_uv)
        force = np.sum([[0.0,0.0], 2*magnitude*np.array(sum_uv)],axis=0)
        force_dict[outer_arc[i]] = np.add(force_dict[outer_arc[i]], force) 

        # i+1 node goes in u-dir
        force = np.sum([[0.0,0.0], magnitude*np.array(udir)],axis=0)
        force_dict[outer_arc[i+1]] = np.add(force_dict[outer_arc[i+1]], force) 

    # update location of node 
    pos = nx.get_node_attributes(G,'pos')
    
    for node in force_dict:
        G.node[node]['pos'] = d_pos(pos[node],force_dict[node],dt)

# Save nx Graph in pickled form for plotting later
    
    if t % 1 == 0: 
        file_name = 't' + str(round(t)) 
        nx.write_gpickle(G,file_name + '.pickle')
        np.save(file_name,circum_sorted)
