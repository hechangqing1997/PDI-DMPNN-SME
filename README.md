#先验数据引入的定向消息传递网络+结构单元掩蔽解释（PDI-DMPNN+SME）联用算法
#prior data introducing directed message passing neural network 
and substructure mask explanation (PDI-DMPNN+SME) coupling algorithm

本项目基于图神经网络DMPNN改进开发，感谢原作者团队。
本项目已授权中国软件著作权（2025SR0349859）

————————————————————————————————————
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
