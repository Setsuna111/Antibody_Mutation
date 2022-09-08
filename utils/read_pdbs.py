import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from protein_attribute import aa_codes, BACKBONE_ATOMs
import Bio
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser


def read_pdb_3D(pdb_file_path, chain_index_list=None, is_backbone=False):
    store_all = []
    with open(pdb_file_path) as protein:
        for lines in protein:
            if 'ATOM' in lines:
                # todo: if is_backbone
                lines = lines.split()
                store_all.append(list(map(float, lines[6:9])))
    return store_all


def show_pdb(store_all, chain_index_list=None):
    x, y, z = zip(*store_all)
    fig = plt.Figure()
    ax = Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(ax)
    ax.plot(x, y, z)
    plt.show()


def read_pdb_seq(pdb_file_path):
    seq = ''
    seq_index = 0
    with open(pdb_file_path) as protein:
        for lines in protein:
            if 'ATOM' in lines:
                lines = lines.split()
                aa_index = lines[5]
                if aa_index == seq_index:
                    continue
                seq_index = aa_index
                seq = seq + aa_codes[lines[3]]  # 乐 这tm是读入原子序列
    return seq


def read_pdb_Biopython(pdb_file_path):
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    structure = parser.get_structure("1CBW", pdb_file_path)
    # print(data)
    model = structure[0]
    print(model.child_dict)
    chain_F = model["F"]
    # chains = list(model.get_chains())


if __name__ == '__main__':
    pdb_file_path = '../data/SKEMPIv2/PDBs_filtered/1CBW.pdb'
    store_all = read_pdb_3D(pdb_file_path)
    print(read_pdb_seq(pdb_file_path))
    # show_pdb(store_all)

    read_pdb_Biopython(pdb_file_path)
