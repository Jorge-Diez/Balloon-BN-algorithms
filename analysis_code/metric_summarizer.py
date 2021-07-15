import os
import numpy
from natsort import natsorted
from distutils.dir_util import copy_tree
from scipy import stats
import pandas as pd
import numpy as np
import statistics
import pickle

NUM_EXPERIMENTS = 88


#Where BOLD data will be copied
METRIC_SAVE_LOCATION = r"D:\My Documents\MUIA\TFM\CasualFMRI\METRICS\GOBLIN"




def progressbar(acc, total, total_bar, exp, run):

    
    
    frac = acc/total
    filled_progbar = round(frac*total_bar)
    print('\r', '#'*filled_progbar + '-'*(total_bar-filled_progbar) + str(exp) + "  " +  str(run) + "|| Exp count: " + str(acc), '[{:>7.2%}]'.format(frac), end='')



  
  
exp_count = 0


for (directory) in natsorted(os.listdir(METRIC_SAVE_LOCATION)):

    main_experiment_directory = os.path.join(METRIC_SAVE_LOCATION, directory)
    
    for (subdirectory) in natsorted(os.listdir(main_experiment_directory)):
        
        experiment_directory = os.path.join(main_experiment_directory, subdirectory)
        all_runs_metrics = []
        
        
        
        for (run) in natsorted(os.listdir(experiment_directory)):
            
            dic_metrics = os.path.join(experiment_directory, run) + "\\mean_metrics.pickle"

            with open( dic_metrics , 'rb') as handle:
                run_mean_metrics = (pickle.load(handle))

        
            all_runs_metrics.append(run_mean_metrics)
        
        
            
        all_connection_tp = []
        all_connection_tn = []
        all_connection_fp = []
        all_connection_fn = []
        all_c_sensitivity = []
        all_c_specificity = []
        all_c_accuracy = []
        all_direction_tp = []
        all_direction_tn = []
        all_direction_fp = []
        all_direction_fn = []
        all_d_sensitivity = []
        all_d_specificity = []
        all_d_accuracy = []
        all_d_smith = []
        
        
            
        for runmetric in all_runs_metrics:
            all_connection_tp.append(runmetric["connection_tp"])
            all_connection_tn.append(runmetric["connection_tn"])
            all_connection_fp.append(runmetric["connection_fp"])
            all_connection_fn.append(runmetric["connection_fn"])
            all_c_sensitivity.append(runmetric["c-sensitivity"])
            all_c_specificity.append(runmetric["c-specificity"])
            all_c_accuracy.append(runmetric["c-accuracy"])
            all_direction_tp.append(runmetric["direction_tp"])
            all_direction_tn.append(runmetric["direction_tn"])
            all_direction_fp.append(runmetric["direction_fp"])
            all_direction_fn.append(runmetric["direction_fn"])
            all_d_sensitivity.append(runmetric["d-sensitivity"])
            all_d_specificity.append(runmetric["d-specificity"])
            all_d_accuracy.append(runmetric["d-accuracy"])
            all_d_smith.append(runmetric["d-smith"])
            
        
        mean_connection_tp = sum(all_connection_tp) / len(all_connection_tp)
        mean_connection_tn = sum(all_connection_tn) / len(all_connection_tn)
        mean_connection_fp = sum(all_connection_fp) / len(all_connection_fp)
        mean_connection_fn = sum(all_connection_fn) / len(all_connection_fn)
        mean_c_sensitivity = sum(all_c_sensitivity) / len(all_c_sensitivity)
        mean_c_specificity = sum(all_c_specificity) / len(all_c_specificity)
        mean_c_accuracy = sum(all_c_accuracy) / len(all_c_accuracy)
        mean_direction_tp = sum(all_direction_tp) / len(all_direction_tp)
        mean_direction_tn = sum(all_direction_tn) / len(all_direction_tn)
        mean_direction_fp = sum(all_direction_fp) / len(all_direction_fp)
        mean_direction_fn = sum(all_direction_fn) / len(all_direction_fn)
        mean_d_sensitivity = sum(all_d_sensitivity) / len(all_d_sensitivity)
        mean_d_specificity = sum(all_d_specificity) / len(all_d_specificity)
        mean_d_accuracy = sum(all_d_accuracy) / len(all_d_accuracy)
        mean_d_smith = sum(all_d_smith) / len(all_d_smith)

        txt_results = experiment_directory + "\\results.txt"

        with open(txt_results, 'x') as f:
            f.write("experiment: {}\n".format(subdirectory))
            f.write("mean_connection_tp : {}\n".format(mean_connection_tp))
            f.write("mean_connection_tn : {}\n".format(mean_connection_tn))
            f.write("mean_connection_fp : {}\n".format(mean_connection_fp))
            f.write("mean_connection_fn : {}\n".format(mean_connection_fn))
            f.write("mean_c_sensitivity : {}\n".format(mean_c_sensitivity))
            f.write("mean_c_specificity : {}\n".format(mean_c_specificity))
            f.write("mean_c_accuracy : {}\n".format(mean_c_accuracy))
            f.write("mean_direction_tp : {}\n".format(mean_direction_tp))
            f.write("mean_direction_tn : {}\n".format(mean_direction_tn))
            f.write("mean_direction_fp : {}\n".format(mean_direction_fp))
            f.write("mean_direction_fn : {}\n".format(mean_direction_fn))
            f.write("mean_d_sensitivity : {}\n".format(mean_d_sensitivity))
            f.write("mean_d_specificity : {}\n".format(mean_d_specificity))
            f.write("mean_d_accuracy : {}\n".format(mean_d_accuracy))
            f.write("mean_d_smith : {}\n".format(mean_d_smith))            
            f.close()
            
        progressbar(exp_count, NUM_EXPERIMENTS, 50, subdirectory, run)
        exp_count += 1

