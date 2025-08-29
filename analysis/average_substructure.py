# -*- codeing = utf-8 -*-
# @Time : 2024/6/13 20:30
# @Author : hcq
# @File : average_substructure.py
# @Software : PyCharm

import pandas as pd

# 读取CSV文件
input_file = 'new_brics.csv'  # 请将文件名替换为你的文件名
df = pd.read_csv(input_file)

# 对substructure_smiles相同的行进行归类，并计算normalized_attribution_data的平均值
grouped_df = df.groupby('substructure_smiles')['normalized_attribution_data'].mean().reset_index()

# 生成新的CSV文件
output_file = 'average_'+input_file  # 生成的文件名
grouped_df.to_csv(output_file, index=False)

print(f"新的CSV文件已生成：{output_file}")
