import numpy as np
import pandas as pd
import os
from torch.utils.data.dataset import Dataset
import torch
from read_pdbs import read_pdb_3D, read_pdb_seq


class SKEMPIV2Dataset(Dataset):
    def __init__(self, data_df, is_train):
        super(SKEMPIV2Dataset, self).__init__()
        self.data_df = data_df
        self.data_batches = len(data_df['PDB_id'])
        self.is_train = is_train

    def __len__(self):
        return self.data_batches

    def rand(self, a=0, b=1):
        return np.random.rand() * (b - a) + a

    def __getitem__(self, index):
        index = index % self.data_batches
        sample_info = self.data_df.iloc[index].values
        PDB_id, chain1, chain2, mutate_info, _, ddG = sample_info
        PDB_file_path = '../data/SKEMPIv2/PDBs_filtered/' + PDB_id + '.pdb'
        # 此处可以选择读入原始坐标，可以在此处进行特征提取，也可以传出后进行特征提取，也可以直接读入外部文件提取好的特征，此处涉及到是否传出突变信息
        # 注意传出后是对数据进行矩阵运算，这里还涉及到是否使用collate，所以特征提取要小心
        complex_full_3D = read_pdb_3D(PDB_file_path, is_backbone=False)
        complex_bb_3D = read_pdb_3D(PDB_file_path, is_backbone=True)

        return PDB_file_path, ddG, chain1


def skempi_dataset_collate(batch):
    pass
    # images_a = []
    # images_b = []
    # bboxes_a = []
    # bboxes_b = []
    # for image_a, image_b, target_a, target_b in batch:
    #     images_a.append(image_a)
    #     bboxes_a.append(target_a)
    #     images_b.append(image_b)
    #     bboxes_b.append(target_b)
    # images_a = np.array(images_a)
    # bboxes_a = np.array(bboxes_a)
    # images_b = np.array(images_b)
    # bboxes_b = np.array(bboxes_b)
    # return images_a, images_b, bboxes_a, bboxes_b
