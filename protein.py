#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 12:45:54 2020

@author: smylonas
"""

import os, numpy as np
import pybel
from utils import simplify_dms


class Protein:
    def __init__(self, prot_file, f):
        prot_id, prot_ext = prot_file.split('/')[-1].split('.')
        if prot_ext != 'pdb':
	    raise IOError('Protein file should be .pdb')
        
        self.mol = next(pybel.readfile(prot_ext,prot_file)) 

        surfpoints_file = os.path.join(prot_file.rsplit('/',1)[0],prot_id+'.surfpoints')
        os.system('dms '+prot_file+' -d 0.2 -n -o '+surfpoints_file)
        if not os.path.exists(surfpoints_file):
            raise Exception('probably DMS not installed')
        
        self.surf_points, self.surf_normals = simplify_dms(surfpoints_file,f)  
        os.remove(surfpoints_file)
            
        self.atom_coords = np.array([atom.coords for atom in self.mol.atoms])

        self.binding_sites = []
        
        with open(prot_file,'r') as f:    
            lines = f.readlines()
        self.atom_lines = [line for line in lines if line[:4]=='ATOM']            
        
        
    def _surfpoints_to_atoms(self,surfpoints):
        coords = self.atom_coords
        close_atoms = np.zeros(len(surfpoints),dtype=int)
        for p,surf_coord in enumerate(surfpoints):
            dist = np.sqrt(np.sum((coords-surf_coord)**2,axis=1))
            close_atoms[p] = np.argmin(dist)
        
        return np.unique(close_atoms)
        
    def add_bsite(self,cluster):   # cluster -> tuple: (surf_points,scores)
        bsite_atom_idxs = self._surfpoints_to_atoms(cluster[0])

        bsite = {
            'coords': self.atom_coords[bsite_atom_idxs],
            'score': np.average(cluster[1]),
            'atom_ids': [int(self.atom_lines[idx].split()[1]) for idx in bsite_atom_idxs]
        }
         
        self.binding_sites.append(bsite)
        
    def sort_bsites(self):
        avg_scores = np.array([bsite['score'] for bsite in self.binding_sites])
        sorted_idxs = np.flip(np.argsort(avg_scores),axis=0)
        self.binding_sites = [self.binding_sites[idx] for idx in sorted_idxs]
        
    
