# DTL-DephosphoSite - Deep Transfer Learning based approach to predict dephosphorylation sites

DTL-DephosphoSite is a transfer learning-based approach that employs deep learning to predict the dephosphorylation sites of S,T, and Y.

DTL-DephosphoSite is a deep learning-based method for Dephosphorylation sites of S,T, and Y. It utilizes phosphorylation data through a transfer learning method to address prediction and present scarce dephosphorylation data.It is implemented using Keras (version 2.2.4) and Tensorflow (version 1.15) backend and has been tested on both Windows and Linux OS.
 

## Requirements

The script requires the following dependencies:

- **Python**: 3.9.18  
- **Libraries**:
  - matplotlib: 3.9.2
  - Biopython: 1.78
  - keras: 2.10.0
  - numpy: 1.26.4
  - scikit-learn: 1.5.1
  - pandas: 2.2.2
  - imbalanced-learn: (optional, for training?)

## Installation

To install these dependencies, use the following command:

```bash
pip install matplotlib==3.9.2 biopython==1.78 keras==2.10.0 numpy==1.26.4 scikit-learn==1.5.1 pandas==2.2.2
```
Alternatively, if you are using Conda, you can create an environment with these versions:

```bash
conda create --name phosphatase_analysis python=3.9 matplotlib=3.9.2 biopython=1.78 keras=2.10.0 numpy=1.26.4 scikit-learn=1.5.1 pandas=2.2.2
conda activate phosphatase_analysis
```

## Running on CPU or GPU
To run in CPU, installation of Tensorflow and Keras will suffice. However, to run in GPU, further Tensorflow-gpu and keras-gpu must be installed. 
Tensorflow GPU and Keras GPU version utilizes cuda cores in our GPU for faster training time. However, running in GPU is not mandatory.
 
## Dataset
The dataset is in FASTA format. Both training and testing datasets are provided, and they are independent (one does not include others). The training dataset for positive and negative are X.fasta and train_X.fasta, respectively. The testing dataset for positive and negative are test_X.fasta and test_X.fasta, respectively. The training dataset is made available so that future models can be trained for comparison purposes.
 
## Model
The model learned through transfer learning from the phosphorylation data for residues ST and Y is provided. The ComDephos_ST.h5 and ComDephos_Y are the optimized models for ST and Y, respectively.

## Code
Independent Testing code is provided. The model provided can be used to predict input dephosphorylation of S, T, and Y sites for the given window sequences. The code will take input of window size 31.

# Prediction for given test dataset (Procedure)
  - Download test datasets, test_Pos_ST.fasta and test_Neg_ST.fasta(from dataset folder), and python code test_model.py.
    Keep them in the same folder as model files.
  - Run test_model.py, and you will get the output mentioned in our research.
  - In Linux code will be $python3 test_model.py
  
## Prediction for your dataset
If you want to use DTL-DephosphoSite to predict dephosphorylation sites in the protein of your interest, prepare your dataset in the same format as the test dataset in FASTA format. 
This model works for window **size 33 only**, meaning you should provide **16 residues downstream and 16 residues upstream** for the residue of your interest.

The general format for your dataset should be:

>sp|Q4KWH8|PLCH1_HUMAN%730%755

PKKQLILKVISGQQLPKPPDSMFGDSGEIIDPFVEVEIIGLPVDCCKDQTR
