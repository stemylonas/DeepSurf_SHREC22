#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 13:39:49 2019

@author: smylonas
"""

import numpy as np
from scipy.spatial.distance import euclidean
from sklearn.cluster import KMeans


mapping = {
    '0.6000': ('H','H'),
    '1.1000': ('HA','H'),
    '1.3590': ('HE1','H'),
    '1.3870': ('HA','H'),
    '1.4090': ('HD1','H'),
    '1.4590': ('HE1','H'),
    '1.4870': ('HB1','H'),
    '1.6612': ('O','O'),
    '1.7210': ('OG','O'),
    '1.8240': ('N','N'),
    '1.9080': ('C','C'),
    '2.0000': ('S','S'),
    }


def pqr_to_pdb(input_file,out_file):

    with open(input_file,'r') as f:
        lines = f.readlines()
    
    pdb_lines = []
    for line in lines:
        parts = line.strip().split()
        atom_id = parts[1]
        x = parts[5]
        y = parts[6]
        z = parts[7]
        radius = parts[9]
        
        if radius == '0.0000':
            continue
        atom_type, element = mapping[radius]
        
        new_line = 'ATOM' + (7-len(atom_id))*' ' + atom_id + '  '
        new_line += atom_type + (4-len(atom_type))*' ' + 'PRO     1    '
        new_line += (8-len(x))*' ' + x +  (8-len(y))*' ' + y + (8-len(z))*' ' + z + '  1.00  0.00'
        new_line +=  11*' ' + element + '\n'
        
        pdb_lines.append(new_line)
    
    with open(out_file,'w') as f:
        f.writelines(pdb_lines)
    
    
def save_bsites(pockets_found, input_pqr_file, out_file):

    with open(input_pqr_file,'r') as f:
        input_pqr_lines = f.readlines()

    for i,pocket in enumerate(pockets_found):
        pocket_id = i+1
        id_str = str(pocket_id)+'.0000' if pocket_id < 10 else str(pocket_id)+'.000' 

        for i,line in enumerate(input_pqr_lines):
            if int(line.split()[1]) in pocket['atom_ids']:
                input_pqr_lines[i] = line[:59]+id_str+line[65:]
        
    with open(out_file,'w') as f:
        f.writelines(input_pqr_lines)
    

def readSurfPoints(surf_file):
    with open(surf_file,'r') as f:
        lines = f.readlines()
    
    lines = [l for l in lines if len(l.split())>7]
    if len(lines)>100000:
        print('Large size')
        return
    if len(lines)==0:
        print('Empty file')
        return
    
    coords = np.zeros((len(lines),3))
    normals = np.zeros((len(lines),3))
    for i,l in enumerate(lines):
        parts = l.split()
        try:
            coords[i,0] = float(parts[3])
            coords[i,1] = float(parts[4])
            coords[i,2] = float(parts[5])
            normals[i,0] = float(parts[8])
            normals[i,1] = float(parts[9])
            normals[i,2] = float(parts[10])
        except:
            coords[i,0] = float(parts[2][-8:])
            coords[i,1] = float(parts[3])
            coords[i,2] = float(parts[4])
            normals[i,0] = float(parts[7])
            normals[i,1] = float(parts[8])
            normals[i,2] = float(parts[9])
            
    return coords, normals


def simplify_dms(init_surf_file, factor):
    
    coords, normals = readSurfPoints(init_surf_file)
    
    if factor == 1:
        return coords, normals

    nPoints = len(coords)
    nCl = nPoints//factor
    
    kmeans = KMeans(n_clusters=nCl,max_iter=300,n_init=1).fit(coords)
    point_labels = kmeans.labels_
    centers = kmeans.cluster_centers_
    cluster_idx,freq = np.unique(point_labels,return_counts=True)
    if len(cluster_idx)!=nCl:  # need to be removed
        print('error')

    idxs = []
    for cl in cluster_idx:
        cluster_points_idxs = np.where(point_labels==cl)[0]
        closest_idx_to_center = np.argmin([euclidean(centers[cl],coords[idx]) for idx in cluster_points_idxs])
        idxs.append(cluster_points_idxs[closest_idx_to_center])
    
    return coords[idxs], normals[idxs]


def rotation(n):
    if n[0]==0.0 and n[1]==0.0:
        if n[2]==1.0:
            return np.identity(3)
        elif n[2]==-1.0:
            Q = np.identity(3)
            Q[0,0] = -1
            return Q
        else:
            print('not possible')
        
    rx = -n[1]/np.sqrt(n[0]*n[0]+n[1]*n[1])
    ry = n[0]/np.sqrt(n[0]*n[0]+n[1]*n[1])
    rz = 0
    th = np.arccos(n[2])
    
    q0 = np.cos(th/2)
    q1 = np.sin(th/2)*rx
    q2 = np.sin(th/2)*ry
    q3 = np.sin(th/2)*rz
               
    Q = np.zeros((3,3))
    Q[0,0] = q0*q0+q1*q1-q2*q2-q3*q3
    Q[0,1] = 2*(q1*q2-q0*q3)
    Q[0,2] = 2*(q1*q3+q0*q2)
    Q[1,0] = 2*(q1*q2+q0*q3)
    Q[1,1] = q0*q0-q1*q1+q2*q2-q3*q3
    Q[1,2] = 2*(q3*q2-q0*q1)
    Q[2,0] = 2*(q1*q3-q0*q2)
    Q[2,1] = 2*(q3*q2+q0*q1)
    Q[2,2] = q0*q0-q1*q1-q2*q2+q3*q3
     
    return Q                                                                              

 
