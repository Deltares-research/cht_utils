# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 16:30:44 2022

@author: ormondt
"""

import numpy as np
from hydrolib.core.io.net.models import Network
from pathlib import Path

class SnapWaveMesh():
    
    def __init__(self):
        
        self.network = Network()

    def to_file(self, file_name):
        p = Path(file_name)
        self.network.to_file(p)

    def index_to_file(self, file_name):
        # Save index_quadtree_in_snapwave file
        # Add 1 because of Fortran in which indexing starts at 1
        # This file is needed when coupling to SFINCS
        file = open(file_name, "wb")
        file.write(np.int32(self.nr_quadtree_cells))
        file.write(np.int32(self.index_quadtree_in_snapwave + 1))
        file.close()
                
    def from_quadtree(self, qtr, mask,
                      file_name=None,
                      index_file_name=None):

        inode  = -1        
        nfaces = -1
        face_nodes = np.zeros((4, 4*qtr.nr_cells), dtype=int) - 1
        index0_quadtree_in_snapwave = np.zeros(qtr.nr_cells, dtype=int) - 1
        
        self.nr_quadtree_cells = qtr.nr_cells
        
        # Determine faces        
        for ip in range(qtr.nr_cells):

            if mask[ip] > 0:

                inode += 1 
                index0_quadtree_in_snapwave[inode] = ip
                
                mu  = qtr.mu[ip]
                mu1 = qtr.mu1[ip]
                mu2 = qtr.mu2[ip]
                nu  = qtr.nu[ip]
                nu1 = qtr.nu1[ip]
                nu2 = qtr.nu2[ip]
                
                # Set inactive in case of mask == 0
                if mu1 >= 0:
                    if mask[mu1] == 0:
                        mu1 = -1
                if mu2 >= 0:
                    if mask[mu2] == 0:
                        mu2 = -1
                if nu1 >= 0:
                    if mask[nu1]==0:
                        nu1 = -1
                if nu2 >= 0:
                    if mask[nu2] == 0:
                        nu2 = -1
                
                mnu  = 0
                mnu1 = -1
                
                # Find neighbors (above-right)
                
                # Try going the right
                if mu == 0:
                    # Same level right
                    if mu1 >= 0:
                        if qtr.nu[mu1] == 0:
                            # Same level above right
                            if qtr.nu1[mu1] >= 0:
                                # and it exists
                                if mask[qtr.nu1[mu1]]>0:
                                    mnu  = 0
                                    mnu1 = qtr.nu1[mu1]
                        elif qtr.nu[mu1] == 1:   
                            # Finer above-right
                            if qtr.nu1[mu1] >= 0:
                                # And it exists
                                if mask[qtr.nu1[mu1]]>0:
                                    mnu  = 1
                                    mnu1 = qtr.nu1[mu1]
                        elif qtr.nu[mu1] == -1:   
                            # Coarser above-right
                            if qtr.nu1[mu1] >= 0:
                                # And it exists
                                if mask[qtr.nu1[mu1]]>0:
                                    mnu  = -1
                                    mnu1 = qtr.nu1[mu1]
                elif mu == -1:
                    # Coarser to the right
                    if mu1 >= 0:
                        if qtr.nu[mu1] == 0:
                            # Same level above right
                            if qtr.nu1[mu1] >= 0:
                                # And it exists
                                if mask[qtr.nu1[mu1]]>0:
                                    mnu  = -1
                                    mnu1 = qtr.nu1[mu1]
                        elif qtr.nu[mu1] == 1:
                            # Finer above right
                            if qtr.nu1[mu1] >= 0:
                                # And it exists
                                if mask[qtr.nu1[mu1]]>0:
                                    mnu  = 0
                                    mnu1 = qtr.nu1[mu1]
                        elif qtr.nu[mu1] == -1:
                            # Coarser above right
                            if qtr.nu1[mu1] >= 0:
                                # And it exists
                                if mask[qtr.nu1[mu1]]>0:
                                    mnu  = -2
                                    mnu1 = qtr.nu1[mu1]
                else:
                    # Finer to the right
                    if mu2 >= 0:
                        if qtr.nu[mu2] == 0:
                            # Same level above right
                            if qtr.nu1[mu2] >= 0:
                                # And it exists
                                mnu  = 1
                                mnu1 = qtr.nu1[mu2]
                        elif qtr.nu[mu2] == 1:
                            # Finer above right
                            if qtr.nu1[mu2] >= 0:
                                # And it exists
                                mnu  = 2
                                mnu1 = qtr.nu1[mu2]
                        else:
                            # Finer above right
                            if qtr.nu1[mu2] >= 0:
                                # And it exists
                                mnu  = 0
                                mnu1 = qtr.nu1[mu2]
                
                # Okay, found all the neighbors!
                    
                # Now let's see what sort of cell this is
                if mu==0 and nu==0 and mnu==0 and mu1>=0 and nu1>=0 and mnu1>=0:
                    # Type 1
                    # Most normal cell possible
                    nfaces += 1
                    face_nodes[0, nfaces] = ip
                    face_nodes[1, nfaces] = mu1
                    face_nodes[2, nfaces] = mnu1
                    face_nodes[3, nfaces] = nu1
                elif (mu==1 and nu==0 and mnu==0):
                    # Type 2
                    if mu1>=0 and mu2>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu1
                        face_nodes[2, nfaces] = mu2
                    if mu2>=0 and mnu1>=0 and nu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = mu2
                        face_nodes[1, nfaces] = mnu1
                        face_nodes[2, nfaces] = nu1
                    if mu2>=0 and nu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu2
                        face_nodes[2, nfaces] = nu1
                elif (mu==0 and nu==0 and mnu==1):
                    # Type 3
                    if mu1>=0 and mnu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu1
                        face_nodes[2, nfaces] = mnu1
                    if nu1>=0 and mnu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mnu1
                        face_nodes[2, nfaces] = nu1
                elif (mu==0 and nu==1 and mnu==0):
                    # Type 4
                    if mu1>=0 and nu2>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu1
                        face_nodes[2, nfaces] = nu2
                    if mu1>=0 and mnu1>=0 and nu2>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = mu1
                        face_nodes[1, nfaces] = mnu1
                        face_nodes[2, nfaces] = nu2
                    if nu1>=0 and nu2>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = nu1
                        face_nodes[2, nfaces] = nu2
                elif (mu==1 and nu==0 and mnu==1):
                    # Type 5
                    if mu1>=0 and mu2>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu1
                        face_nodes[2, nfaces] = mu2
                    if mu2>=0 and mnu1>=0 and nu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu2
                        face_nodes[2, nfaces] = mnu1
                        face_nodes[3, nfaces] = nu1
                elif (mu==0 and nu==1 and mnu==1):
                    # Type 6
                    if mu1>=0 and mnu1>=0 and nu2>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu1
                        face_nodes[2, nfaces] = mnu1
                        face_nodes[3, nfaces] = nu2
                    if nu1>=0 and nu2>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = nu2
                        face_nodes[2, nfaces] = nu1
                elif (mu==1 and nu==1 and (mnu==1 or mnu==0)):
                    # Type 7 and 8
                    if mu1>=0 and mu2>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu1
                        face_nodes[2, nfaces] = mu2
                    if mu2>=0 and nu2>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu2
                        face_nodes[2, nfaces] = nu2
                    if mu2>=0 and nu2>=0 and mnu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = mu2
                        face_nodes[1, nfaces] = mnu1
                        face_nodes[2, nfaces] = nu2
                    if nu1>=0 and nu2>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = nu2
                        face_nodes[2, nfaces] = nu1
                elif (mu==-1 and nu==0 and not odd(qtr.n[ip])):
                    # Type 9
                    if mu1>=0 and nu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu1
                        face_nodes[2, nfaces] = nu1
        #         elif (mu==-1 and nu==0 and mnu==0 and not odd(qtr.n[ip]))
        #             # Type 9
        #             if mu1>0 and nu1>0
        #                 nfaces += 1
        #                 face_nodes[0, nfaces] = ip
        #                 face_nodes[1, nfaces] = mu1
        #                 face_nodes[2, nfaces] = nu1
        #             end
                elif (mu==-1 and nu==-1 and mnu==-1):
                    # Type 10
                    if mu1>=0 and mnu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu1
                        face_nodes[2, nfaces] = mnu1
                    if mnu1>=0 and nu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mnu1
                        face_nodes[2, nfaces] = nu1
                elif (mu==-1 and nu==-1 and mnu==0):
                    # Type 11
                    if mu1>=0 and mnu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu1
                        face_nodes[2, nfaces] = mnu1
                    if mnu1>=0 and nu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mnu1
                        face_nodes[2, nfaces] = nu1
                elif mu==0 and nu==-1 and mnu==-1 and odd(qtr.m[ip]):
                    # Type 12
                    if mu1>=0 and mnu1>=0 and nu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu1
                        face_nodes[2, nfaces] = mnu1
                        face_nodes[3, nfaces] = nu1
                elif (mu==-1 and nu==0 and mnu==-1 and odd(qtr.n[ip])):
                    # Type 13
                    if mu1>=0 and mnu1>=0 and nu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu1
                        face_nodes[2, nfaces] = mnu1
                        face_nodes[3, nfaces] = nu1
                elif (mu==0 and nu==-1 and not odd(qtr.m[ip])):
                    # Type 14
                    if mu1>=0 and nu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu1
                        face_nodes[2, nfaces] = nu1
                elif (mu==1 and nu==-1 and mnu==0):
                    # Type 15
                    if mu1>=0 and mu2>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu1
                        face_nodes[2, nfaces] = mu2
                    if mu2>=0 and mnu1>=0 and nu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = mu2
                        face_nodes[1, nfaces] = mnu1
                        face_nodes[2, nfaces] = nu1
                    if mu2>=0 and nu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu2
                        face_nodes[2, nfaces] = nu1
                elif (mu==-1 and nu==-1 and mnu==-2):
                    # Type 16
                    if mu1>=0 and mnu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu1
                        face_nodes[2, nfaces] = mnu1
                    if mnu1>=0 and nu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mnu1
                        face_nodes[2, nfaces] = nu1
                elif (mu==0 and nu==0 and mnu==-1):
                    # Type 16
                    if mu1>=0 and nu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu1
                        face_nodes[2, nfaces] = nu1
                    if mu1>=0 and nu1>=0 and mnu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = mu1
                        face_nodes[1, nfaces] = mnu1
                        face_nodes[2, nfaces] = nu1
                elif (mu==0 and nu==-1 and mnu==0):
                    # Type 17
                    if mu1>=0 and mnu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu1
                        face_nodes[2, nfaces] = mnu1
                    if mnu1>=0 and nu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mnu1
                        face_nodes[2, nfaces] = nu1
                elif (mu==1 and nu==1 and mnu==2):
                    # Type 17
                    if mu1>=0 and mu2>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu1
                        face_nodes[2, nfaces] = mu2
                    if mu2>=0 and mnu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu2
                        face_nodes[2, nfaces] = mnu1
                    if nu2>=0 and mnu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mnu1
                        face_nodes[2, nfaces] = nu2
                    if nu1>=0 and nu2>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = nu1
                        face_nodes[2, nfaces] = nu2
                elif (mu==1 and nu==1 and mnu==2):
                    # Type 18
                    if mu1>=0 and mu2>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu1
                        face_nodes[2, nfaces] = mu2
                    if mu2>=0 and mnu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu2
                        face_nodes[2, nfaces] = mnu1
                    if nu2>=0 and mnu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mnu1
                        face_nodes[2, nfaces] = nu2
                    if nu1>=0 and nu2>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = nu1
                        face_nodes[2, nfaces] = nu2
                elif (mu==-1 and nu==1 and mnu==0):
                    # Type 19
                    if mu1>=0 and nu2>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu1
                        face_nodes[2, nfaces] = nu2
                    if nu2>=0 and mnu1>=0 and mu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = mu1
                        face_nodes[1, nfaces] = mnu1
                        face_nodes[2, nfaces] = nu2
                    if nu1>=0 and nu2>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = nu2
                        face_nodes[2, nfaces] = nu1
                elif (mu==-1 and nu==0 and mnu==0 and odd(qtr.n[ip])):
                    # Type 20
                    if mu1>=0 and mnu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mu1
                        face_nodes[2, nfaces] = mnu1
                    if mnu1>=0 and nu1>=0:
                        nfaces += 1
                        face_nodes[0, nfaces] = ip
                        face_nodes[1, nfaces] = mnu1
                        face_nodes[2, nfaces] = nu1
        
        nfaces += 1
        nnodes = inode + 1

        # Remove empty faces
        face_nodes = face_nodes[:, 0:nfaces]
        
        # Remove trailing indices
        index0_quadtree_in_snapwave = index0_quadtree_in_snapwave[0:nnodes]      

        # There may be nodes left that are not connected to any faces
        # We need to get rid of these
        self.index_quadtree_in_snapwave = np.empty(nnodes, dtype=int)
        index_snapwave_in_quadtree = np.empty(qtr.nr_cells, dtype=int)
        j = 0
        for inode in range(nnodes):
            # Check if this node is part of a face
            # Index of the quadtree cell isw
            isw = index0_quadtree_in_snapwave[inode]
            if np.size(np.where(face_nodes==isw)[0])>0:
                self.index_quadtree_in_snapwave[j] = isw
                index_snapwave_in_quadtree[isw] = j
                j += 1
        nnodes = j        
        self.index_quadtree_in_snapwave = self.index_quadtree_in_snapwave[0:j]             
        
        node_x  = np.zeros(nnodes)
        node_y  = np.zeros(nnodes)
        node_z  = np.zeros(nnodes)
        face_x  = np.zeros(nfaces)
        face_y  = np.zeros(nfaces)
        
        for inode in range(nnodes):
            isw = self.index_quadtree_in_snapwave[inode]
            node_x[inode] = qtr.x[isw]
            node_y[inode] = qtr.y[isw]
            node_z[inode] = qtr.z[isw]
        
        # Indices in faces are still for full quadtree
        # Change to indices of reduced mesh
        for iface in range(nfaces):
            for k in range(4):
                if face_nodes[k, iface]>=0:
                    face_nodes[k, iface] = index_snapwave_in_quadtree[face_nodes[k, iface]]
                else:
                    # This is a triangle
                    pass

        # Set face_x and face_y to 0.0. These are not used anyway.                
        for iface in range(nfaces):
            face_x[iface] = 0.0
            face_y[iface] = 0.0
                        
        # Edges
        edge_nodes = np.zeros((2, nfaces*4), dtype=int)
        iedge = 0
        # First add all edges for each face
        for iface in range(nfaces):
            edge_nodes[0, iedge] = face_nodes[0, iface]
            edge_nodes[1, iedge] = face_nodes[1, iface]
            iedge += 1
            edge_nodes[0, iedge] = face_nodes[1, iface]
            edge_nodes[1, iedge] = face_nodes[2, iface]
            iedge += 1
            if face_nodes[3, iface] == -1:
                # Triangle
                edge_nodes[0, iedge] = face_nodes[2, iface]
                edge_nodes[1, iedge] = face_nodes[0, iface]
                iedge += 1
            else:
                # Rectangle
                edge_nodes[0, iedge] = face_nodes[2, iface]
                edge_nodes[1, iedge] = face_nodes[3, iface]
                iedge += 1
                edge_nodes[0, iedge] = face_nodes[3, iface]
                edge_nodes[1, iedge] = face_nodes[0, iface]
                iedge += 1
        edge_nodes = edge_nodes[:, 0:iedge]
        # Now sort and get rid of double edges 
        edge_nodes = np.sort(edge_nodes, axis=0)        
        edge_nodes = np.sort(edge_nodes, axis=1)        
        edge_nodes = np.unique(edge_nodes, axis=1)        
                        
        # And now we make the ugrid mesh
        face_nodes[np.where(face_nodes<0)]     = -1000 # one will be added
        self.network._mesh2d.mesh2d_node_x     = node_x
        self.network._mesh2d.mesh2d_node_y     = node_y
        self.network._mesh2d.mesh2d_node_z     = node_z
        self.network._mesh2d.mesh2d_edge_nodes = np.transpose(edge_nodes)
        self.network._mesh2d.mesh2d_face_nodes = np.transpose(face_nodes)
        self.network._mesh2d.mesh2d_face_x     = face_x
        self.network._mesh2d.mesh2d_face_y     = face_y
        
        if file_name:
            self.to_file(file_name)
            
        if index_file_name:
            self.index_to_file(index_file_name)
            

    def plot(self, ax, edgecolor="k", facecolor="r", fill=True):
        
        itriangle = np.where(self.network._mesh2d.mesh2d_face_nodes[:, 3]<0)
        self.network._mesh2d.mesh2d_face_nodes[itriangle, 3] = self.network._mesh2d.mesh2d_face_nodes[itriangle, 0]

        import matplotlib.patches as patches

        xface = self.network._mesh2d.mesh2d_node_x[self.network._mesh2d.mesh2d_face_nodes]
        yface = self.network._mesh2d.mesh2d_node_y[self.network._mesh2d.mesh2d_face_nodes]
        nfaces = np.shape(self.network._mesh2d.mesh2d_face_nodes)[0]
        for i in range(nfaces):
            path = [ [xface[i, 0], yface[i, 0]],
                     [xface[i, 1], yface[i, 1]],
                     [xface[i, 2], yface[i, 2]],
                     [xface[i, 3], yface[i, 3]] ] 
            p = patches.Polygon(path,
                                facecolor=facecolor,
                                edgecolor=edgecolor,
                                fill=fill)
            ax.add_patch(p)        

def odd(num):
    if (num % 2) == 1:  
        return True
    else:  
        return False
