#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 16:04:47 2020

@author: smylonas
"""

from network import Network
import argparse, os, time
from protein import Protein
from bsite_extraction import Bsite_extractor
from tqdm import tqdm
from utils import pqr_to_pdb, save_bsites


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--input_path', '-inp', required=True, help='directory with the input protein files')
    parser.add_argument('--model_path', '-mp', required=True, help='directory of models')
    parser.add_argument('--model', '-m', choices=['orig','lds'], default='orig', help='select model')
    parser.add_argument('--output', '-o', required=True, help='name of the output directory')
    parser.add_argument('--f', type=int, default=10, help='parameter for the simplification of points mesh')
    parser.add_argument('--T', type=float, default=0.9, help='ligandability threshold')
    parser.add_argument('--batch', type=int, default=32, help='batch size')
    parser.add_argument('--voxel_size', type=float, default=1.0, help='size of voxel in angstrom')

    return parser.parse_args()


args = parse_args()

if not os.path.exists(args.input_path):
    raise IOError('%s does not exist.' % args.input_path)
if not os.path.exists(args.model_path):
    raise IOError('%s does not exist.' % args.model_path)
if not os.path.exists(args.output):
    os.makedirs(args.output)

st = time.time()

nn = Network(args.model_path,args.model,args.voxel_size)

extractor = Bsite_extractor(args.T)

protein_ids = os.listdir(args.input_path)

for protein_id in tqdm(protein_ids):
    pqr_file = os.path.join(args.input_path,protein_id,'structure.pqr')
    temp_pdb_file = os.path.join(args.input_path,protein_id+'.pdb')
    pqr_to_pdb(pqr_file,temp_pdb_file)

    protein = Protein(temp_pdb_file, args.f)
    os.remove(temp_pdb_file)
    
    lig_scores = nn.get_lig_scores(protein,args.batch)
    
    extractor.extract_bsites(protein,lig_scores)
    
    save_bsites(protein.binding_sites, pqr_file, os.path.join(args.output,protein_id+'.pqr'))
    
print('Total time:',time.time()-st)
