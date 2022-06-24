#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 16:16:56 2020

@author: smylonas
"""

from network import Network
import argparse, os
from protein import Protein
from bsite_extraction import Bsite_extractor
from utils import pqr_to_pdb, save_bsites

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--input_file', '-p', required=True, help='input file (pqr)')
    parser.add_argument('--model_path', '-mp', required=True, help='directory of models')
    parser.add_argument('--model', '-m', choices=['orig','lds'], default='orig', help='select model')
    parser.add_argument('--output', '-o', required=True, help='name of the output directory')
    parser.add_argument('--f', type=int, default=10, help='parameter for the simplification of points mesh')
    parser.add_argument('--T', type=float, default=0.9, help='ligandability threshold')
    parser.add_argument('--batch', type=int, default=32, help='batch size')
    parser.add_argument('--voxel_size', type=float, default=1.0, help='size of voxel in angstrom')

    return parser.parse_args()


args = parse_args()

if not os.path.exists(args.input_file):
    raise IOError('%s does not exist.' % args.input_file)
if not os.path.exists(args.model_path):
    raise IOError('%s does not exist.' % args.model_path)
if not os.path.exists(args.output):
    os.makedirs(args.output)

prot_id, prot_ext = args.input_file.split('/')[-1].split('.')
if prot_ext != 'pqr':
    raise IOError('Input file should be .pqr')

temp_pdb_file = os.path.join(args.input_path.rsplit('/',1)[0], prot_id+'.pdb')
pqr_to_pdb(args.input_file,temp_pdb_file)

prot = Protein(temp_pdb_file, args.f)
os.remove(temp_pdb_file)

nn = Network(args.model_path,args.model,args.voxel_size)

lig_scores = nn.get_lig_scores(prot,args.batch)

extractor = Bsite_extractor(args.T)

extractor.extract_bsites(prot,lig_scores)

save_bsites(prot.binding_sites, args.input_file, os.path.join(args.output,prot_id+'.pqr'))
    
