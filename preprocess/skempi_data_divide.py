import pandas as pd
import numpy as np
import os
from shutil import copyfile
import sklearn
from sklearn.model_selection import train_test_split

if __name__ == '__main__':
    data = pd.read_csv('../data/SKEMPIv2/skempi_v2_filtered.csv')
    PDB_id_list = data['PDB_id']
    key_list = ['PDB_id', 'chain1', 'chain2', 'Mutation_info', 'Binding_site_location', 'ddG']
    data_train = pd.DataFrame(columns=key_list)
    data_val = pd.DataFrame(columns=key_list)
    data_test = pd.DataFrame(columns=key_list)
    region_list = ['INT', 'SUR', 'RIM', 'COR', 'SUP']

    for region in region_list:
        data_region = data[data['Binding_site_location'] == region]

        data_region_x, data_region_test = train_test_split(data_region, test_size=0.1, train_size=0.9
                                                           , shuffle=True, random_state=301)
        data_region_train, data_region_val = train_test_split(data_region_x, test_size=2/9, train_size=7/9
                                                              , shuffle=True, random_state=301)
        data_train = data_train.append(data_region_train)
        data_val = data_val.append(data_region_val)
        data_test = data_test.append(data_region_test)

    # 再把多位点的按比例分进去
    data_region = data[(data['Binding_site_location'] != 'INT') & (data['Binding_site_location'] != 'SUR') &
                       (data['Binding_site_location'] != 'RIM') & (data['Binding_site_location'] != 'COR') &
                       (data['Binding_site_location'] != 'SUP')]
    data_region_x, data_region_test = train_test_split(data_region, test_size=0.1, train_size=0.9
                                                       , shuffle=True, random_state=301)
    data_region_train, data_region_val = train_test_split(data_region_x, test_size=2 / 9, train_size=7 / 9
                                                          , shuffle=True, random_state=301)
    data_train = data_train.append(data_region_train)
    data_val = data_val.append(data_region_val)
    data_test = data_test.append(data_region_test)
    print(len(data), len(data_train), len(data_val), len(data_test))

    data_train.to_csv('../data/SKEMPIv2/skempi_v2_train.csv', index=None)
    data_val.to_csv('../data/SKEMPIv2/skempi_v2_val.csv', index=None)
    data_test.to_csv('../data/SKEMPIv2/skempi_v2_test.csv', index=None)


