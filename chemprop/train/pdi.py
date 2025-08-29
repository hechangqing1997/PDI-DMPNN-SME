# -*- codeing = utf-8 -*-
# @Time : 2023/12/2 13:42
# @Author : hcq
# @File : PDI.py
# @Software : PyCharm

import csv
from chemprop.args import CommonArgs


def csv_to_dict(csv_filename):
    with open(csv_filename, 'r') as file:
        reader = csv.reader(file)
        num_columns = get_num_columns(csv_filename)
        dict_list = [{} for _ in range(num_columns - 1)]

        for i in range(1, num_columns):
            dict_name = dict_list[i - 1]
            file.seek(0)  # 重新将迭代指针移动到文件开头
            for row in reader:
                if row:
                    smiles = row[0]
                    value = row[i]
                    dict_name[smiles] = value
    return dict_list

def batch_prior_data(batch_smiles):
    solvent_filename=CommonArgs.solvent_filename
    solute_filename=CommonArgs.solute_filename
    
    batch_prior_data=[]
    solvent_dict_list = csv_to_dict(solvent_filename)
    solute_dict_list = csv_to_dict(solute_filename)

    for sublist in batch_smiles:
        solvent_smile=sublist[0]
        solute_smile= sublist[1]

        solvent_prior_data=search_prior_data(solvent_dict_list, solvent_smile)
        solute_prior_data = search_prior_data(solute_dict_list, solute_smile)

        difference_prior_data=compute_absolute_difference(solvent_prior_data, solute_prior_data)
        batch_prior_data.append(difference_prior_data)

    return batch_prior_data

def compute_absolute_difference(solvent_prior_data, solute_prior_data):
    if len(solvent_prior_data) != len(solute_prior_data):
        raise ValueError("溶剂和溶质先验列表长度不一致")

    difference = [abs(x - y) for x, y in zip(solvent_prior_data, solute_prior_data)]
    return difference


def get_num_columns(csv_filename):
    with open(csv_filename, 'r') as file:
        csv_reader = csv.reader(file)
        # 获取文件中第一行的列数
        first_row = next(csv_reader)
        num_columns = len(first_row)
    return num_columns




def search_prior_data(dict_list, search_smiles):
    try:
        data_list = []
        for dict_name in dict_list:
            if search_smiles in dict_name:
                data_list.append(dict_name[search_smiles])

        data_list = [float(item) for item in data_list]
        
        return data_list

    except Exception as e:
        print(f"发生了一个错误：{e}")
    return []

