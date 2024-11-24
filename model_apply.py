#%%
# # -*- coding: utf-8 -*-
"""
Change res in line with "res = 'Y'", to get the chosen residue results
"""

# ---------------------------------------------------------------------------- #
# Model Evaluation

import matplotlib.pyplot as plt
# plt.ion()  # Enable interactive mode
import matplotlib
matplotlib.use('TkAgg')
from Bio import SeqIO
import keras
import numpy as np
from numpy import array
from sklearn.metrics import roc_curve, auc, classification_report
from keras.models import load_model
import pandas as pd

# import Bio
# import sklearn
# print("matplotlib:", matplotlib.__version__)
# print("Biopython:", Bio.__version__)
# print("keras:", keras.__version__)
# print("numpy:", np.__version__)
# print("scikit-learn:", sklearn.__version__)
# print("pandas:", pd.__version__)

# ---------------------------------------------------------------------------- #

# Configuration

# Input parameters
res = 'ST'  # Specify residue type: 'Y', 'S/T'
condition = 'arsenite_decreased'  # Specify condition: treatment, 'increased', 'decreased'
input_fasta = f"dataset/{condition}_output_{res}.fasta"  # Input FASTA file
model_path = f"ComDephos_{res}.h5"  # Path to model
output_csv = f"model_output/{condition}_{res}.csv"  # Output file for predictions

# Model input window parameters
if res == 'Y': # Update to match model input - Y: 27, ST: 31
    win = 27 
elif res == 'ST':
    win = 31
win_size = 33
cut_off = int((33 - win) / 2)

# Alphabet for sequence encoding
alphabet = 'ARNDCQEGHILKMFPSTWYV*'
char_to_int = dict((c, i) for i, c in enumerate(alphabet))

# ---------------------------------------------------------------------------- #
# Helper functions

def encode_sequence(seq):
    """
    Encode the sequence as integers.
    """
    seq = seq[cut_off:-cut_off]
    for char in seq:
        if char not in alphabet:
            return None
    return [char_to_int[char] for char in seq]

# ---------------------------------------------------------------------------- #
# Process FASTA file

# Prepare data for prediction
x_test = []
seq_names = []

for seq_record in SeqIO.parse(input_fasta, "fasta"):
    encoded = encode_sequence(seq_record.seq)
    if encoded:
        x_test.append(encoded)
        seq_names.append(seq_record.id)  # Save sequence names for mapping results

x_test = array(x_test)  # Convert to NumPy array

# ---------------------------------------------------------------------------- #
# Load model and predict

model = load_model(model_path)
predictions = model.predict(x_test)
predicted_probabilities = predictions[:, 1]  # Extract probabilities for dephosphorylation

# ---------------------------------------------------------------------------- #
# Save results

results = pd.DataFrame({
    "Sequence ID": seq_names,
    "Dephosphorylation Probability": predicted_probabilities
})

threshold = 0.5  # Define threshold
results["Dephosphorylation"] = results["Dephosphorylation Probability"].apply(
    lambda x: "yes" if x >= threshold else "no"
)

results.to_csv(output_csv, index=False)
print(f"Predictions classified and saved to {output_csv}")
