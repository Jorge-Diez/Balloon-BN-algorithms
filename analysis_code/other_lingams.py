import os
import numpy
from natsort import natsorted
from distutils.dir_util import copy_tree
from scipy import stats
import pandas as pd
import numpy as np
import lingam
from lingam.utils import make_dot

def matrix_transform_triplets (data_matrix):
    
    tripleterinos = []
    
    for node_from,column in enumerate(data_matrix.T):
        for node_to, connection_value in enumerate(column):
        
            if (connection_value != 0):
                print("from node: {}  to node:  {}  there is a connection strength of: {}".format(node_from, node_to, connection_value))
                triplet = []
                triplet_from = "C" + str(node_from + 1)
                triplet.append(triplet_from)
                triplet.append('-->')
                triplet_to = "C" + str(node_to + 1)
                triplet.append(triplet_to)
                
                tripleterinos.append(triplet)
                
    
    return tripleterinos






colnames=["C1", "C2", "C3", "C4", "C5"] 
data_stringloc = r"D:\My Documents\MUIA\TFM\CasualFMRI\NETSIM_Balloon_Data\Resting_State\Resting-State_5min\run1\BOLD_DATA\FILTERED\subject1bold_filtered.csv"
data = r"D:\My Documents\MUIA\TFM\CasualFMRI\NETSIM_Balloon_Data\Resting_State\Resting-State_5min\run1\BOLD_DATA\FILTERED\subject26bold_filtered.csv"
data = pd.read_csv (data, names = colnames, header=None)



model = lingam.DirectLiNGAM()
model.fit(data)

adj_matrix = model.adjacency_matrix_

print(adj_matrix)
print("   \n") 


newtriplets = matrix_transform_triplets(adj_matrix)








