 # -*- codeing = utf-8 -*-
# @Time : 2024/6/3 13:35
# @Author : hcq
# @File : attribution_calculate.py
# @Software : PyCharm
import csv
import numpy as np
from rdkit import Chem

csv_filename='with+without_solvent-brics_mol_graph_group_for_brics.csv'
npy_filename='solvent-brics_mol_graph_smask_for_brics.npy'
sme_mask_column='mol solvent'
with_sme_column = 'with_sme'
without_sme_column = 'without_sme'
new_csv_filename = 'new_brics.csv'

def get_sub_molecule(mol, atom_indices):
    mol=Chem.MolFromSmiles(mol)
    # 创建一个空的分子对象
    sub_mol = Chem.RWMol()
    # 创建一个映射字典来保存新原子索引
    atom_map = {}

    # 将选定的原子添加到新分子中
    for idx in atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        new_idx = sub_mol.AddAtom(atom)
        atom_map[idx] = new_idx

    # 将键添加到新分子中
    for bond in mol.GetBonds():
        start_atom = bond.GetBeginAtomIdx()
        end_atom = bond.GetEndAtomIdx()
        if start_atom in atom_indices and end_atom in atom_indices:
            sub_mol.AddBond(atom_map[start_atom], atom_map[end_atom], bond.GetBondType())

    # 转换为分子对象
    sub_mol = sub_mol.GetMol()
    sub_smiles = Chem.MolToSmiles(sub_mol)
    return sub_smiles


# Open original CSV file for reading
with open(csv_filename, mode='r') as file:
    csv_reader = csv.DictReader(file)

    # Open new CSV file for writing
    with open(new_csv_filename, mode='w', newline='') as new_file:
        fieldnames = csv_reader.fieldnames + ['attribution_data', 'normalized_attribution_data','substructure_smiles']
        csv_writer = csv.DictWriter(new_file, fieldnames=fieldnames)
        csv_writer.writeheader()

        # Process each row in original CSV
        for line_num, row in enumerate(csv_reader, start=1):
            with_sme_data = np.float32(row.get(with_sme_column, 0))
            without_sme_data = np.float32(row.get(without_sme_column, 0))
            attribution_data = without_sme_data - with_sme_data
            normalized_attribution_data = (np.exp(attribution_data) - np.exp(-attribution_data)) / (np.exp(attribution_data) + np.exp(-attribution_data))

            sme_mask_smiles=row.get(sme_mask_column, 0)
            sme_mask_substructures = np.load(npy_filename, allow_pickle=True)
            atom_indices=sme_mask_substructures[line_num-1]
            substructure_smiles=get_sub_molecule(sme_mask_smiles, atom_indices)

            # Add calculated data to the row and write it to the new CSV file
            row['attribution_data'] = attribution_data
            row['normalized_attribution_data'] = normalized_attribution_data
            row['substructure_smiles']=substructure_smiles
            csv_writer.writerow(row)

            print('第'+str(line_num)+'行已经处理完毕！')

print("已完成所有的分子片段处理！")