import os
import numpy
from natsort import natsorted
from distutils.dir_util import copy_tree
from scipy import stats
import pandas as pd
import numpy as np
import statistics
import pickle
import rpy2.robjects as ro
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from pygobnilp.gobnilp import Gobnilp
import contextlib



def obtain_all_connections (data_direction, nnodes) :
    all_lines = []
    with open (data_direction, "r") as file:
        all_lines = file.readlines()
        
    first_index = nnodes + 1
    last_index = all_lines.index("stop\n")
    
    clean_lines = [ data_line[0:len(data_line)-1] for data_line in all_lines[first_index:last_index]  ]
    
    return clean_lines



def triplet_creator_all(line_data) :
    all_tripleterinos = []
    
    for line in line_data:
        tripleto = []
        only_important = line[5:10]
        node_split = only_important.rsplit(" ")
        tripleto.append(node_split[0])
        tripleto.append("-->")
        tripleto.append(node_split[1])
        all_tripleterinos.append(tripleto)
        
    
    return all_tripleterinos



#PARA HACER EL GAMBVA MIENTRAS PRUEBAS CON OTROS METODOS


colnames=["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"] 
data = r"D:\My Documents\MUIA\TFM\CasualFMRI\NETSIM_Balloon_Data_Nocycle_deletablegoblin\Boxcar_Empirical_Resting_State\Boxcar_Empirical_Resting-State_excellentconditions\run8\BOLD_DATA\FILTERED\subject46bold_filtered.csv"
data = pd.read_csv (data, names = colnames, header=None)

nnodes = 8



m = Gobnilp()    
returnedshit = m.learn(data, data_type='continuous',score='BGe',palim=None, plot=False, output_stem = "stemerino", output_ext = ('plain',))


all_lines = obtain_all_connections('stemerino.plain', nnodes)
all_triplets = triplet_creator_all(all_lines)