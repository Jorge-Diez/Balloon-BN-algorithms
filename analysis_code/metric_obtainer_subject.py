import os
import numpy
from natsort import natsorted
from distutils.dir_util import copy_tree
from scipy import stats
import pandas as pd
import numpy as np
import statistics
import pickle

NUM_EXPERIMENTS = 52800



"""
I am fully aware that some parts are clearly not fully efficient, I am just running on my last fumes 
and do not have the time to program everything efficiently
"""

#where is all the bold data
SIM_DATA_FOLDER = r"D:\My Documents\MUIA\TFM\CasualFMRI\NETSIM_ALGORITHM_RESULTS\LINGAM"

#Where BOLD data will be copied
METRIC_SAVE_LOCATION = r"D:\My Documents\MUIA\TFM\CasualFMRI\METRICS\LINGAM"




smith = [["C1", "-->", "C2"], ["C1", "-->", "C5"], ["C2", "-->", "C3"], ["C3", "-->", "C4"], ["C4", "-->", "C5"]]
smithcycle = [["C1", "-->", "C2"], ["C2", "-->", "C1"], ["C1", "-->", "C5"], ["C2", "-->", "C3"], ["C3", "-->", "C4"], ["C4", "-->", "C5"]]

smithstim = [["C4", "-->", "C3"], ["C4", "-->", "C5"], ["C3", "-->", "C2"], ["C5", "-->", "C1"], ["C2", "-->", "C1"]]
smithstimcycle = [["C4", "-->", "C3"], ["C4", "-->", "C5"], ["C3", "-->", "C2"], ["C2", "-->", "C3"], ["C5", "-->", "C1"], ["C2", "-->", "C1"]]

experimental = [["C1", "-->", "C3"], ["C1", "-->", "C5"], ["C2", "-->", "C3"], ["C2", "-->", "C5"], ["C3", "-->", "C4"], ["C3", "-->", "C6"], ["C3", "-->", "C8"], ["C6", "-->", "C7"]]
experimentalstim = [["C1", "-->", "C3"], ["C1", "-->", "C5"], ["C2", "-->", "C3"], ["C2", "-->", "C5"], ["C3", "-->", "C4"], ["C3", "-->", "C6"], ["C3", "-->", "C8"], ["C6", "-->", "C7"]]
experimentalcycle = [["C1", "-->", "C3"], ["C1", "-->", "C5"], ["C5", "-->", "C1"], ["C2", "-->", "C3"], ["C2", "-->", "C5"], ["C3", "-->", "C4"], ["C4", "-->", "C3"], ["C3", "-->", "C6"], ["C3", "-->", "C8"], ["C6", "-->", "C7"]]
experimentalstimcycle = [["C1", "-->", "C3"], ["C1", "-->", "C5"], ["C5", "-->", "C1"], ["C2", "-->", "C3"], ["C2", "-->", "C5"], ["C3", "-->", "C4"], ["C4", "-->", "C3"], ["C3", "-->", "C6"], ["C3", "-->", "C8"], ["C6", "-->", "C7"]]

smith_nr_nodes = 5
experimental_nr_nodes = 8





def progressbar(acc, total, total_bar, exp, run):

    
    
    frac = acc/total
    filled_progbar = round(frac*total_bar)
    print('\r', '#'*filled_progbar + '-'*(total_bar-filled_progbar) + str(exp) + "  " +  str(run) + "|| Subject count: " + str(acc), '[{:>7.2%}]'.format(frac), end='')








def separate_strings(line_data) :
    
    quitar_numero = line_data.rsplit(".")[1]
    raw_text = quitar_numero[1:len(quitar_numero)]
    
    return raw_text.rsplit(" ")



def obtain_all_connections (data_direction) :
    all_lines = []
    with open (data_direction, "r") as file:
        all_lines = file.readlines()
        
    first_index = all_lines.index("Graph Edges:\n") + 1
    last_index = all_lines.index("Graph Attributes:\n") - 1
    
    clean_lines = [ data_line[0:len(data_line)-1] for data_line in all_lines[first_index:last_index]  ]
    
    return clean_lines


def detect_experiment (folder):
    
    if (folder == "Resting_State"):
        return smith, smith_nr_nodes
    elif (folder == "Cycle_Resting_State"):
        return smithcycle, smith_nr_nodes
    elif (folder == "Boxcar_Resting_State"):
        return smithstim, smith_nr_nodes
    elif (folder == "Cycle_Boxcar_Resting_State"):
        return smithstimcycle, smith_nr_nodes
    elif (folder == "Empirical_Resting_State"):
        return experimental, experimental_nr_nodes
    elif (folder == "Boxcar_Empirical_Resting_State"):
        return experimentalstim, experimental_nr_nodes
    elif (folder == "Cycle_Empirical_Resting_State"):
        return experimentalcycle, experimental_nr_nodes
    elif (folder == "Cycle_Boxcar_Empirical_Resting_State"):
        return experimentalstimcycle, experimental_nr_nodes
    
    
    
def check_connection_positive (datatriplet, groundtruth):
    """
    Done here because it makes code less jumbled

    """
    isconnection = 0
    

    for truth_triplet in groundtruth:
        if (  (datatriplet[0] == truth_triplet[0] and datatriplet[2] == truth_triplet[2]) or
              (datatriplet[0] == truth_triplet[2] and datatriplet[2] == truth_triplet[0]) ):
            isconnection = 1
            break
        
    return isconnection




def check_direction (datatriplet, groundtruth):
    """
    only for one single direction

    """
    isconnection = 0
    

    for truth_triplet in groundtruth:
        if (  (datatriplet[0] == truth_triplet[0] and datatriplet[2] == truth_triplet[2]) ):
            isconnection = 1
            break
        
    return isconnection


def create_metric (data, truth, metrics, nnodes):
    
    nr_connections = len(truth)
    tp_positions = [] #positions where there are tp in data
    
    c_tp_count = 0
    c_fp_count = 0
    
    for i,triplet_d in enumerate(data):
        
        if (check_connection_positive(triplet_d, truth)):
            c_tp_count += 1
            tp_positions.append(i)
        else:
            c_fp_count += 1
     
    metrics['connection_tp'].append(c_tp_count)
    metrics['connection_fp'].append(c_fp_count)
        
    #How many false negatives? FN mean that they are NOT in data, but they ARE in truth
    false_negatives = nr_connections - (metrics['connection_tp'][-1])
    metrics['connection_fn'].append(false_negatives)
    
    #true negatives means that connections are NOT in data and NOT in truth as well
    #so, the total count of true negatives is the total amount of negatives - false positives?
    possible_connections = (nnodes*nnodes) - nnodes #matrix - self connections
    all_true_negatives = possible_connections - len(truth)
    metrics['connection_tn'].append(all_true_negatives - metrics['connection_fp'][-1])
    
    
    metrics['c-sensitivity'].append(metrics['connection_tp'][-1] / nr_connections)
    metrics['c-specificity'].append(metrics['connection_tn'][-1] / (possible_connections - nr_connections))
    metrics['c-accuracy'].append((metrics['connection_tp'][-1] + metrics['connection_tn'][-1]) / ((all_true_negatives) + nr_connections))
    
    
    #d-smith
    smith_hits = 0
    d_tp_count = 0
    d_fp_count = 0
    
    
    for j,direction_triplet in enumerate(data):
        
        
        if direction_triplet[1] == "-->":
            if (check_direction(direction_triplet, truth)):
                
                d_tp_count += 1
                
                if (j in tp_positions): # for the smith direction metric
                    smith_hits += 1
                
                
            else:
                d_fp_count += 1
                
    metrics['direction_tp'].append(d_tp_count)
    metrics['direction_fp'].append(d_fp_count)
     
    #remember, we are still using nr connectiond and such because all of our connections in truth have direction
        
    false_negatives_direction = nr_connections - metrics['direction_tp'][-1]
    metrics['direction_fn'].append(false_negatives_direction)
        
    possible_connections = (nnodes*nnodes) - nnodes #matrix - self connections
    all_true_negatives = possible_connections - len(truth)
    metrics['direction_tn'].append(all_true_negatives - metrics['direction_fp'][-1] )
    
    
    metrics['d-sensitivity'].append(metrics['direction_tp'][-1] / nr_connections)
    metrics['d-specificity'].append(metrics['direction_tn'][-1] / (possible_connections - nr_connections))
    metrics['d-accuracy'].append((metrics['direction_tp'][-1] + metrics['direction_tn'][-1]) / ((all_true_negatives) + nr_connections))
    
    if (smith_hits == 0):
        metrics['d-smith'].append(0)
    else:
        metrics['d-smith'].append(smith_hits / len(tp_positions))
  
    
    
    
    
    
    
exp_count = 0


for (directory) in natsorted(os.listdir(SIM_DATA_FOLDER)):

    main_experiment_directory = os.path.join(SIM_DATA_FOLDER, directory)
    
    for (subdirectory) in natsorted(os.listdir(main_experiment_directory)):
        
        experiment_directory = os.path.join(main_experiment_directory, subdirectory)
        
        for (run) in natsorted(os.listdir(experiment_directory)):
            
            experiment_directory_run = os.path.join(experiment_directory, run)
            ground_truth, nr_nodes = detect_experiment(directory)

            all_metrics = { # tp: true positive, tn: true negative, likewise for f = false
                "data_origin": experiment_directory_run,
                "connection_tp": [], 
                "connection_tn": [],
                "connection_fp": [],
                "connection_fn": [],
                "c-sensitivity": [],
                "c-specificity": [],
                "c-accuracy": [],
                "direction_tp": [], 
                "direction_tn": [],
                "direction_fp": [],
                "direction_fn": [],
                "d-sensitivity": [],
                "d-specificity": [],
                "d-accuracy": [],
                "d-smith": []
                #"cycle_tp": [], 
                #"cycle_tn": [],
                #"cycle_fp": [],
                #"cycle_fn": [],
                #"cyc-sensitivity": [],
                #"cyc-specificity": [],
                #"cyc-accuracy": [],
                }
            
            mean_metrics = {
                "connection_tp": [], 
                "connection_tn": [],
                "connection_fp": [],
                "connection_fn": [],
                "c-sensitivity": [],
                "c-specificity": [],
                "c-accuracy": [],
                "direction_tp": [], 
                "direction_tn": [],
                "direction_fp": [],
                "direction_fn": [],
                "d-sensitivity": [],
                "d-specificity": [],
                "d-accuracy": [],
                "d-smith": []
                }
            
            for subject_result in natsorted(os.listdir(experiment_directory_run)):
                
                result_txt = os.path.join(experiment_directory_run, subject_result)
                connections = obtain_all_connections(result_txt)
                sub_out = subject_result.rsplit(".txt")[0]
                subnum = int(sub_out.rsplit("subject")[1])
                
                all_triplets = []
                for clean_line in connections:
                    all_triplets.append(separate_strings(clean_line))
                    
                
                create_metric(all_triplets, ground_truth, all_metrics, nr_nodes)
                    
                progressbar(exp_count, NUM_EXPERIMENTS, 50, subdirectory, run)
                exp_count += 1
                
                
            
                
            direction_metrics = METRIC_SAVE_LOCATION + "\\" + directory +  "\\" + subdirectory + "\\" + run
            os.makedirs(direction_metrics)
            
            with open( direction_metrics + "\\raw_metrics" + '.pickle', 'wb') as handle:
                pickle.dump(all_metrics, handle, protocol=pickle.HIGHEST_PROTOCOL) 
                
            for key in mean_metrics: 
                data = all_metrics[key]
                mean_metrics[key] = sum(data) / len(data)
                
            with open( direction_metrics + "\\mean_metrics" + '.pickle', 'wb') as handle:
                pickle.dump(mean_metrics, handle, protocol=pickle.HIGHEST_PROTOCOL)
                
            
            
