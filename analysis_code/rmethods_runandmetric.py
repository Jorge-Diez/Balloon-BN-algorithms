# Analyze all of the data in a folder with the GES algorithm
import os
import numpy
from natsort import natsorted
import pandas as pd
import rpy2.robjects as ro
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import pickle

NUM_EXPERIMENTS = 52800

#where is all the bold data
DATA_FOLDER = r"D:\My Documents\MUIA\TFM\CasualFMRI\NETSIM_Balloon_Data"



METRIC_SAVE_LOCATION_MMHC = r"D:\My Documents\MUIA\TFM\CasualFMRI\METRICS\MMHC"
METRIC_SAVE_LOCATION_TABU = r"D:\My Documents\MUIA\TFM\CasualFMRI\METRICS\TABU"

r = robjects.r
r['source']('mmhc_rcode.R')


funcion_mmhc = robjects.globalenv['mmhc_calc']
funcion_tabu = robjects.globalenv['tabu_calc']


smith = [["V1", "-->", "V2"], ["V1", "-->", "V5"], ["V2", "-->", "V3"], ["V3", "-->", "V4"], ["V4", "-->", "V5"]]
smithcycle = [["V1", "-->", "V2"], ["V2", "-->", "V1"], ["V1", "-->", "V5"], ["V2", "-->", "V3"], ["V3", "-->", "V4"], ["V4", "-->", "V5"]]

smithstim = [["V4", "-->", "V3"], ["V4", "-->", "V5"], ["V3", "-->", "V2"], ["V5", "-->", "V1"], ["V2", "-->", "V1"]]
smithstimcycle = [["V4", "-->", "V3"], ["V4", "-->", "V5"], ["V3", "-->", "V2"], ["V2", "-->", "V3"], ["V5", "-->", "V1"], ["V2", "-->", "V1"]]

experimental = [["V1", "-->", "V3"], ["V1", "-->", "V5"], ["V2", "-->", "V3"], ["V2", "-->", "V5"], ["V3", "-->", "V4"], ["V3", "-->", "V6"], ["V3", "-->", "V8"], ["V6", "-->", "V7"]]
experimentalstim = [["V1", "-->", "V3"], ["V1", "-->", "V5"], ["V2", "-->", "V3"], ["V2", "-->", "V5"], ["V3", "-->", "V4"], ["V3", "-->", "V6"], ["V3", "-->", "V8"], ["V6", "-->", "V7"]]
experimentalcycle = [["V1", "-->", "V3"], ["V1", "-->", "V5"], ["V5", "-->", "V1"], ["V2", "-->", "V3"], ["V2", "-->", "V5"], ["V3", "-->", "V4"], ["V4", "-->", "V3"], ["V3", "-->", "V6"], ["V3", "-->", "V8"], ["V6", "-->", "V7"]]
experimentalstimcycle = [["V1", "-->", "V3"], ["V1", "-->", "V5"], ["V5", "-->", "V1"], ["V2", "-->", "V3"], ["V2", "-->", "V5"], ["V3", "-->", "V4"], ["V4", "-->", "V3"], ["V3", "-->", "V6"], ["V3", "-->", "V8"], ["V6", "-->", "V7"]]

smith_nr_nodes = 5
experimental_nr_nodes = 8




def progressbar(acc, total, total_bar, exp, run):

    
    
    frac = acc/total
    filled_progbar = round(frac*total_bar)
    print('\r', '#'*filled_progbar + '-'*(total_bar-filled_progbar) + str(exp) + "  " +  str(run) + "|| Subject count: " + str(acc), '[{:>7.2%}]'.format(frac), end='')









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
  
    
    








def dat_triplet_transform(pandas_duplet):
    
    tripleterinos = []
      
    for index, row in pandas_duplet.iterrows():
        tripleto = []
        tripleto.append(row['from'])
        tripleto.append('-->')
        tripleto.append(row['to'])
        tripleterinos.append(tripleto)
    
    
    
    
    return tripleterinos































exp_count = 0

for (directory) in natsorted(os.listdir(DATA_FOLDER)):

    main_experiment_directory = os.path.join(DATA_FOLDER, directory)
    
    for (subdirectory) in natsorted(os.listdir(main_experiment_directory)):
        
        experiment_directory = os.path.join(main_experiment_directory, subdirectory)
        
        for (run) in natsorted(os.listdir(experiment_directory)):
            
            experiment_directory_run = os.path.join(experiment_directory, run)
            bold_data_directory = os.path.join(experiment_directory_run, "BOLD_DATA/FILTERED")
            
            
            ground_truth, nr_nodes = detect_experiment(directory)
            
            
            
            all_metrics_mmhc = { # tp: true positive, tn: true negative, likewise for f = false
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
            
            mean_metrics_mmhc = {
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
            
            
            
            
            all_metrics_tabu = { # tp: true positive, tn: true negative, likewise for f = false
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
            
            mean_metrics_tabu = {
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
            
            

            
            
            
            for (subject) in natsorted(os.listdir(bold_data_directory)):
                subject_data = os.path.join(bold_data_directory, subject)
                
                

                result_mmhc = funcion_mmhc(subject_data)
                result_tabu = funcion_tabu(subject_data)

                
                with localconverter(ro.default_converter + pandas2ri.converter):
                    mmhc_data = ro.conversion.rpy2py(result_mmhc)
  
                with localconverter(ro.default_converter + pandas2ri.converter):
                    tabu_data = ro.conversion.rpy2py(result_tabu)
                
                
                all_triplets_mmhc = dat_triplet_transform(mmhc_data)
                all_triplets_tabu = dat_triplet_transform(tabu_data)
                
                create_metric(all_triplets_mmhc, ground_truth, all_metrics_mmhc, nr_nodes)
                create_metric(all_triplets_tabu, ground_truth, all_metrics_tabu, nr_nodes)
                
                exp_count += 1
                progressbar(exp_count, NUM_EXPERIMENTS, 50, subdirectory, run)
                
                
                
            direction_metrics_mmhc = METRIC_SAVE_LOCATION_MMHC + "\\" + directory +  "\\" + subdirectory + "\\" + run
            direction_metrics_tabu = METRIC_SAVE_LOCATION_TABU + "\\" + directory +  "\\" + subdirectory + "\\" + run

            os.makedirs(direction_metrics_mmhc)
            os.makedirs(direction_metrics_tabu)
            
            with open( direction_metrics_mmhc + "\\raw_metrics" + '.pickle', 'wb') as handle:
                pickle.dump(all_metrics_mmhc, handle, protocol=pickle.HIGHEST_PROTOCOL) 
                
            with open( direction_metrics_tabu + "\\raw_metrics" + '.pickle', 'wb') as handle:
                pickle.dump(all_metrics_tabu, handle, protocol=pickle.HIGHEST_PROTOCOL) 
                
                
            for key in mean_metrics_mmhc: 
                data = all_metrics_mmhc[key]
                mean_metrics_mmhc[key] = sum(data) / len(data)
                
            for key in mean_metrics_tabu: 
                data = all_metrics_tabu[key]
                mean_metrics_tabu[key] = sum(data) / len(data)
                
            with open( direction_metrics_mmhc + "\\mean_metrics" + '.pickle', 'wb') as handle:
                pickle.dump(mean_metrics_mmhc, handle, protocol=pickle.HIGHEST_PROTOCOL)
                
            with open( direction_metrics_tabu + "\\mean_metrics" + '.pickle', 'wb') as handle:
                pickle.dump(mean_metrics_tabu, handle, protocol=pickle.HIGHEST_PROTOCOL)
                
                
                
                
                
                
                
                
                
                