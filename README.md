#先验数据引入的定向消息传递网络+结构单元掩蔽解释（PDI-DMPNN+SME）联用算法
#prior data introducing directed message passing neural network 
and substructure mask explanation (PDI-DMPNN+SME) coupling algorithm

# PDI-DMPNN+SME 项目文件架构
## 核心文件结构
chemprop/
├── args.py # 配置参数类定义文件，包含 CommonArgs 等配置类
├── models/ # 模型相关文件夹
│ ├── init.py # 模型组件导出
│ ├── model.py # 主要模型类 MoleculeModel 定义
│ ├── mpn.py # 消息传递网络实现
│ ├── sme_mask.py # SME 掩蔽机制实现
│ └── ffn.py # 前馈神经网络实现
├── train/ # 训练相关文件夹
│ ├── init.py # 训练组件导出
│ └── pdi.py # PDI（先验数据集成）实现
├── data/ # 数据处理相关文件夹
│ └── utils.py # 数据处理工具函数
└── analysis/ # 分析工具文件夹
└── attribution_calculate.py # NDsub计算分析工具
## 主要文件功能说明
### 配置文件
- `args.py`
  - 定义模型所需的各种参数配置
  - 包含 PDI-DMPNN+SME 特有的参数设置
  - 设置数据路径、模型参数等
### 模型核心文件
- `models/model.py`
  - 定义主要模型类 MoleculeModel
  - 实现分子表示的处理和预测
- `models/mpn.py`
  - 实现分子消息传递网络
  - 处理分子图的特征提取
- `models/sme_mask.py`
  - 实现结构单元掩蔽（SME）机制
  - 处理分子结构的掩蔽操作
### 训练相关文件
- `train/pdi.py`
  - 实现先验数据引入（PDI）功能
  - 处理先验知识的整合
### 数据处理文件
- `data/utils.py`
  - 提供数据处理的工具函数
  - 处理分子数据的加载和预处理
### 分析工具
- `analysis/attribution_calculate.py`
  - 计算和分析分子属性
  - 处理 SME 掩蔽前后的数据对比
## 主要功能模块
1. **分子表示学习**：通过 DMPNN 实现
2. **结构掩蔽**：通过 SME 机制实现
3. **先验知识引入**：通过 PDI 机制实现
4. **数据处理和分析**：通过各种工具函数实现

本项目基于图神经网络DMPNN改进开发，感谢原作者团队。
本项目已申请中国软件著作权（2025SR0349859）

——————————————————English——————————————————
### Project File Structure for # PDI-DMPNN+SME
## Core File Structure
```
chemprop/
├── args.py# Configuration parameter class definitions (e.g., CommonArgs)
├── models/# Model implementation directory
│├── __init__.py# Model component exports
│├── model.py# Main model class: MoleculeModel
│├── mpn.py# Message Passing Neural Network (DMPNN)
│├── sme_mask.py# Substructure Masking Enhancement (SME) mechanism
│└── ffn.py# Feed-Forward Network implementation
├── train/# Training framework
│├── __init__.py# Training component exports
│└── pdi.py# Prior Data Integration (PDI) implementation
├── data/# Data processing utilities
│└── utils.py# Molecular data handling functions
└── analysis/# Analytical tools
└── attribution_calculate.py# NDsub attribution analysis tool
```

## Key File Descriptions
### Configuration
- **`args.py`**
- Defines model hyperparameters and configurations
- Includes PDI-DMPNN+SME specific settings (data paths, model architecture, etc.)

### Core Models
- **`models/model.py`**
- Implements `MoleculeModel` class for molecular representation learning and prediction
- **`models/mpn.py`**
- Directed Message Passing Neural Network (DMPNN) for molecular graph feature extraction
- **`models/sme_mask.py`**
- Substructure Masking Enhancement (SME) operations for structural masking

### Training Modules
- **`train/pdi.py`**
- Prior Data Integration (PDI) for knowledge-guided training

### Data Processing
- **`data/utils.py`**
- Molecular data loading, preprocessing, and augmentation utilities

### Analytics
- **`analysis/attribution_calculate.py`**
- Computes molecular attributions (e.g., NDsub) and conducts SME-masking comparative analysis

## Functional Modules
1. **Molecular Representation Learning**
- Implemented via DMPNN in `mpn.py`
2. **Substructure Masking**
- SME mechanism in `sme_mask.py`
3. **Prior Knowledge Integration**
- PDI framework in `pdi.py`
4. **Data Processing & Analysis**
- Utilities in `data/utils.py` and `analysis/attribution_calculate.py`

---
*Note: This project extends the Directed Message Passing Neural Network (DMPNN) with novel SME and PDI modules. Acknowledgment to the original authors.*
**Software Copyright Registered in China (2025SR0349859)**
