from vapory import *
import globals as const
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import itertools
from scipy.spatial import distance

start = 0 
end = 110 

def unit_vector(A,B):
    # Calculate the unit vector from A to B in 3D

    dist = distance.euclidean(A,B)

    if dist < 10e-15:
        dist = 1.0

    return [(B[0]-A[0])/dist,(B[1]-A[1])/dist]

def vector(A,B):
    
    return [(B[0]-A[0]),(B[1]-A[1]), (B[2]-A[2])]

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

def tissue_2d(hex):

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

     
    return G, G, centers, i, circum_sorted, belt

G, H, centers, num_apical_nodes, circum_sorted, belt = tissue_2d(7)

#r = 90 
#theta = np.pi
#xyz = [r*np.cos(theta), r*np.sin(theta), 0]
#
##xyz = [0,0,100]
#
#camera = Camera( 'location', [xyz[0], xyz[1], xyz[2]], 'look_at', [0, 0, 0], 'sky', [0,0,1] )
#light = LightSource( [xyz[0], xyz[1], xyz[2]], 'rgb', [1, 1, 1] )
#background = Background( "rgb", [1,1,1] )
#
for n in range(start,end,10):
    print(n) 
    file_name = 't' + str(n)
   
    G = nx.read_gpickle(file_name + '.pickle')
    circum_sorted = np.load(file_name + '.npy',allow_pickle=True)
   
    pos = nx.get_node_attributes(G,'pos')
    
    mesh_elements = []
    bars = []

    for node in centers:
        nbhrs = circum_sorted[centers.index(node)]
        for i in range(0,len(nbhrs)): 
            # apical triangles
            triangle1 = Triangle(pos[node],pos[nbhrs[i-1]],pos[nbhrs[i]], Texture( Pigment( 'color', [0, 0, 1] ), Finish('phong',0)))
            # apical edges
            cylinder1 = Cylinder(pos[nbhrs[i-1]],pos[nbhrs[i]], 0.25, Texture( Pigment( 'color', [1, 1, 0] )))
            # add to elements for plotting
            mesh_elements.append(triangle1)
            bars.append(cylinder1)

    for i in range(-1,len(const.outer_arc)-1):
            bars.append( Cylinder(pos[const.outer_arc[i-1]],pos[const.outer_arc[i]], 0.5, Texture( Pigment( 'color', [1, 0, 0] ))))
            bars.append( Cylinder(pos[const.outer_arc[i+1]],pos[const.outer_arc[i]], 0.5, Texture( Pigment( 'color', [1, 0, 0] ))))
             
    r = 95 
    theta = -np.pi/2 
    xyz = [r*np.cos(theta), r*np.sin(theta), 0]

    camera = Camera( 'location', [xyz[0], xyz[1], xyz[2]], 'look_at', [0, 0, 0], 'sky', [0,0,1] )
    light = LightSource( [xyz[0], xyz[1], xyz[2]], 'rgb', [1, 1, 1] )
    background = Background( "rgb", [1,1,1] )

    obj_list = [light, background]
    for tri in mesh_elements:
        obj_list.append(tri)
    for bar in bars:
        obj_list.append(bar)

#    scene = Scene( camera, objects = obj_list )
#    scene.render("side%03d.png"%n, width=1000, height=1000)  

    theta = np.pi/2 
    xyz = [np.cos(theta), np.sin(theta), r]
    camera = Camera( 'location', [xyz[0], xyz[0], -r], 'look_at', [0, 0, 0], 'sky', [0,0,1])
    light = LightSource( [xyz[0], xyz[1], -xyz[2]], 'rgb', [1, 1, 1] )
    obj_list[0] = light

    scene = Scene( camera, objects = obj_list )
    scene.render("top%03d.png"%n, width=1000, height=1000)  


