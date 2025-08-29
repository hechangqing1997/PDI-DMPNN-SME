import numpy as np
import pandas as pd
from rdkit import Chem
import csv
import torch
from chemprop.args import CommonArgs,TrainArgs,PredictArgs


sme_mask_npy_filename=CommonArgs.sme_mask_npy_filename
sme_mask_npy = np.load(sme_mask_npy_filename, allow_pickle=True)
test_filename=PredictArgs.test_path

def sme_mask(a_message,batch_number):
    batch_already=TrainArgs.batch_size*(batch_number-1) #已经经历过的分子数
    csv_data = pd.read_csv(test_filename, header=None)   #这个是算上表头的
    molecules = np.array(csv_data.iloc[:, 0])
    molecule_atom_numbers_batch=[]
    
    for k in range (batch_already+1,batch_already+TrainArgs.batch_size+1):
        if k<=sme_mask_npy.shape[0]:
            molecule_atom_number=get_atom_count(molecules[k])
            molecule_atom_numbers_batch.append(molecule_atom_number)

    sme_mask_npy_batch=sme_mask_npy[batch_already:batch_already+TrainArgs.batch_size]
            
    a_message_sme_mask=sme_mask_atoms(a_message,molecule_atom_numbers_batch,sme_mask_npy_batch)   

    return a_message_sme_mask

    
def get_atom_count(smiles: str) -> int:
    # 使用RDKit将SMILES字符串转换为分子对象
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    # 获取分子中的原子数量
    atom_count = mol.GetNumAtoms()
    return atom_count

def sme_mask_atoms(a_message,molecule_atom_numbers_batch,sme_mask_npy_batch):
    num_molecules_batch=sme_mask_npy_batch.shape[0]
    # 第一个部分是单独的(0, 1000)
    parts = [a_message[0].reshape(1, -1)]
    parts_sme=[]
    # 剩余部分根据molecule_atom_numbers_batch进行分割
    start_index = 1
    for atom_num in molecule_atom_numbers_batch:
        end_index = start_index + atom_num
        parts.append(a_message[start_index:end_index])
        start_index = end_index       
    parts_sme.append(parts[0])
    for molecule_num in range(0,num_molecules_batch):            
        rows_to_zero=sme_mask_npy_batch[molecule_num]
        tensor = parts[molecule_num+1]
        for row in rows_to_zero:
            tensor[row, :] = 0   
        parts_sme.append(tensor) 
    
    a_message_sme_mask = torch.cat(parts_sme, dim=0)
    
    return a_message_sme_mask 
    
    
