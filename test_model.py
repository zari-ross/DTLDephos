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
# Residue Selection

res = 'Y'
x_test = []
y_test = []
posit_1 = 1
negat_0 = 0
alphabet = 'ARNDCQEGHILKMFPSTWYV*'
num_classes = 2
win = 27  # Update to match model input
win_size = 33  # Actual window size
cut_off = int((33 - win)/2)

# define a mapping of chars to integers
char_to_int = dict((c, i) for i, c in enumerate(alphabet))
int_to_char = dict((i, c) for i, c in enumerate(alphabet))

# ---------------------------------------------------------------------------- #
# Test DATASET
## For Positive Sequence

def process_positive_sequence():
    data = seq_record.seq
    data = data[cut_off:-cut_off]
    # Integer encode input data
    for char in data:
        if char not in alphabet:
            return
    integer_encoded = [char_to_int[char] for char in data]
    x_test.append(integer_encoded)
    y_test.append(posit_1)

for seq_record in SeqIO.parse("dataset/test_Pos_" + str(res) + ".fasta", "fasta"):
    process_positive_sequence()

# ---------------------------------------------------------------------------- #
## For Negative Sequence

def process_negative_sequence():
    data = seq_record.seq
    data = data[cut_off:-cut_off]
    # Integer encode input data
    for char in data:
        if char not in alphabet:
            return
    integer_encoded = [char_to_int[char] for char in data]
    x_test.append(integer_encoded)
    y_test.append(negat_0)

for seq_record in SeqIO.parse("dataset/test_Neg_" + str(res) + ".fasta", "fasta"):
    process_negative_sequence()

# ---------------------------------------------------------------------------- #
# Changing to array (matrix)

x_test = array(x_test)  # print(x_test.shape) - for debugging from main
test_y1 = array(y_test)
y_test = keras.utils.to_categorical(test_y1, num_classes)

# ---------------------------------------------------------------------------- #
# Load the Model

model_path = "ComDephos_" + str(res) + ".h5"
model = load_model(model_path)
score = model.evaluate(x_test, y_test, verbose=0)  # print(model.input_shape) - for debugging from main
print('Train-val loss:', score[0])
print('Train-val accuracy:', score[1])
acc_train = score[1]

# ---------------------------------------------------------------------------- #
# Scores

from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import confusion_matrix
Y_pred = model.predict(x_test)
t_pred2 = Y_pred[:,1]
Y_pred = (Y_pred > 0.5)
y_pred1 = [np.argmax(y, axis=None, out=None) for y in Y_pred]
y_pred1 = np.array(y_pred1)

print("Matthews Correlation : ",matthews_corrcoef(y_test[:,1],y_pred1))
print("Confusion Matrix : \n",confusion_matrix( y_test[:,1],y_pred1))
mcc_train = matthews_corrcoef(y_test[:,1],y_pred1)

# ---------------------------------------------------------------------------- #
# Sensitivity and Specificity

tn, fp, fn, tp = confusion_matrix(y_test[:,1], y_pred1).ravel()
sp_2_train = tn / (tn + fp)
sn_2_train = tp / (tp + fn)

# ---------------------------------------------------------------------------- #
# ROC

fpr, tpr, _ = roc_curve(y_test[:,1], t_pred2)
roc_auc_train = auc(fpr, tpr)
print("AUC : ", roc_auc_train)
print(classification_report(y_test[:,1], y_pred1))
print("Specificity = ",sp_2_train, " Sensitivity = ",sn_2_train)
plt.figure()
lw = 2
plt.plot(fpr, tpr, color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % roc_auc_train)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC curve for'+str(res))
plt.legend(loc="lower right")
plt.savefig("roc_curve.png")
plt.show()
