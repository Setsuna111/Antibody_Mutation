import pandas as pd
import numpy as np
import os
from shutil import copyfile


def getFlist(path):
    for root, dirs, files in os.walk(file_dir):
        print('root_dir:', root)  # 当前路径
        print('sub_dirs:', dirs)  # 子文件夹
        print('files:', files)    # 文件名称，返回list类型
    return files


if __name__ == '__main__':
    # 处理csv文件，每一行作为一个PDB_sample，并unique处理
    # data = pd.read_csv('data/SKEMPIv2/skempi_v2_raw.csv', sep=';')
    data = pd.read_csv('../data/SKEMPIv2/skempi_v2_changed.csv')

    all_PDB_sample = data['#Pdb']
    unique_PDB_sample = all_PDB_sample.unique()
    print(len(all_PDB_sample))
    print(len(unique_PDB_sample))

    # 筛选PDBs文件，看是否包含了所有sample
    file_dir = '../data/SKEMPIv2/PDBs'
    file_name_list = getFlist(file_dir)

    # count = 0
    for PDB_sample in unique_PDB_sample:
        PDB_id = PDB_sample[:4]
        PDB_file_name = PDB_id + '.pdb'
        if PDB_file_name in file_name_list:
            # count += 1
            continue
        else:
            print(PDB_id, '文件不存在')
    # print(count)  #348个独立的sample 位于345个不同结构，多的那三个是由于重复结构不同链造成的

    # 看看文件类别 发现相同
    pdb_ = 0
    pdb = 0
    mapping_ = 0
    mapping = 0
    for file_name in file_name_list:
        if file_name[:2] == '._':
            if 'mapping' in file_name:
                mapping_ += 1
            else:
                pdb_ += 1
        else:
            if 'mapping' in file_name:
                mapping += 1
            else:
                pdb += 1
    print(pdb, pdb_, mapping, mapping_)

    # 处理csv文件 注意这里没必要分开紧凑表示的突变
    key_list = ['PDB_id', 'chain1', 'chain2', 'Mutation_info', 'Binding_site_location', 'ddG']
    data_new = pd.DataFrame(columns=key_list)
    for i in range(len(all_PDB_sample)):
        PDB_sample = all_PDB_sample[i]
        PDB_id, chain1, chain2 = PDB_sample.split('_')  # 有的链长度过多，需要使用下划线切割
        Mutation_info = data['Mutation(s)_cleaned'].iloc[i]
        Binding_site = data['iMutation_Location(s)'].iloc[i]

        affinity_wt = data['Affinity_wt_parsed'].iloc[i]
        affinity_mut = data['Affinity_mut_parsed'].iloc[i]
        if np.isnan(affinity_wt) or np.isnan(affinity_mut):
            # print('find NAN')
            continue
        dG_wt = (8.314/4184)*(273.15 + 25.0) * np.log(affinity_wt)
        dG_mut = (8.314/4184)*(273.15 + 25.0) * np.log(affinity_mut)
        ddG = dG_mut - dG_wt
        data_new = data_new.append({'PDB_id': str(PDB_id), 'chain1': chain1, 'chain2': chain2, 'Mutation_info': Mutation_info
                           ,'Binding_site_location': Binding_site, 'ddG': ddG}, ignore_index=True)
    print(len(data_new['PDB_id']))
    print(len(data_new['PDB_id'].unique()))
    data_new[data_new.columns[0]] = data_new[data_new.columns[0]].astype('str')
    # data_new.to_csv('../data/SKEMPIv2/skempi_v2_filtered.csv', index=None)

    # 处理一遍pdb文件，放到新的文件夹下 注意判断有无重名 直接用unique的就可以
    pre_file = '../data/SKEMPIv2/PDBs/'
    fin_file = '../data/SKEMPIv2/PDBs_filtered/'
    PDB_id_list = data_new['PDB_id'].unique()
    print(len(PDB_id_list))
    for PDB_id in PDB_id_list:
        pre_pdb = pre_file + PDB_id + '.pdb'
        fin_pdb = fin_file + PDB_id + '.pdb'
        a = os.path.exists(pre_pdb)
        if a == False:
            print('not exist')
        # copyfile(pre_pdb, fin_pdb)

