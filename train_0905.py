import numpy as np
import pandas as pd
import os
from torch.utils.data.dataset import Dataset
from torch.utils.data import DataLoader
import torch
from utils.dataloader import SKEMPIV2Dataset, skempi_dataset_collate
from utils.arguments import get_common_args

if __name__ == '__main__':
    args = get_common_args()
    print(args)

    use_cuda = args.is_cuda
    GPU_indx = 0
    device = torch.device(GPU_indx if use_cuda else "cpu")

    train_path = 'data/SKEMPIv2/skempi_v2_train.csv'
    val_path = 'data/SKEMPIv2/skempi_v2_val.csv'
    test_path = 'data/SKEMPIv2/skempi_v2_test.csv'

    train_df = pd.read_csv(train_path)
    val_df = pd.read_csv(val_path)
    test_df = pd.read_csv(test_path)

    train_dataset = SKEMPIV2Dataset(train_df, is_train=True)
    # gen_train = DataLoader(train_dataset, shuffle=True, batch_size=args.batch_size, pin_memory=False,
    #                        drop_last=True, num_workers=4, collate_fn=skempi_dataset_collate)
    gen_train = DataLoader(train_dataset, shuffle=True, batch_size=args.batch_size, pin_memory=False,
                           drop_last=True, num_workers=4)

    for iteration, batch in enumerate(gen_train):
        print(iteration, batch)




