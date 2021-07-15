import os
import numpy
from natsort import natsorted
from distutils.dir_util import copy_tree
from scipy import stats
import pandas as pd
import numpy as np
import statistics
import pickle
import statistics


NUM_EXPERIMENTS = 308


#Where BOLD data will be copied
METRIC_SAVE_LOCATION = r"D:\My Documents\MUIA\TFM\CasualFMRI\METRICS_STDEV_GOOD"




def progressbar(acc, total, total_bar, exp, run):

    
    
    frac = acc/total
    filled_progbar = round(frac*total_bar)
    print('\r', '#'*filled_progbar + '-'*(total_bar-filled_progbar) + str(exp) + "  " +  str(run) + "|| Exp count: " + str(acc), '[{:>7.2%}]'.format(frac), end='')



  
  
exp_count = 0


for (directory) in natsorted(os.listdir(METRIC_SAVE_LOCATION)):

    main_experiment_directory = os.path.join(METRIC_SAVE_LOCATION, directory)
    
    for (subdirectory) in natsorted(os.listdir(main_experiment_directory)):
        
        experiment_directory = os.path.join(main_experiment_directory, subdirectory)
        
        
        for (concrete_experiment) in natsorted(os.listdir(experiment_directory)):
        
            experiment_directory_full = os.path.join(experiment_directory, concrete_experiment)

            all_runs_metrics = []
        
        
            for (run) in natsorted(os.listdir(experiment_directory_full)):
                
                dic_metrics = os.path.join(experiment_directory_full, run) + "\\mean_metrics.pickle"
    
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
            
        
            stdev_connection_tp = statistics.stdev(all_connection_tp)
            stdev_connection_tn = statistics.stdev(all_connection_tn) 
            stdev_connection_fp = statistics.stdev(all_connection_fp) 
            stdev_connection_fn = statistics.stdev(all_connection_fn) 
            stdev_c_sensitivity = statistics.stdev(all_c_sensitivity) 
            stdev_c_specificity = statistics.stdev(all_c_specificity) 
            stdev_c_accuracy = statistics.stdev(all_c_accuracy) 
            stdev_direction_tp = statistics.stdev(all_direction_tp) 
            stdev_direction_tn = statistics.stdev(all_direction_tn) 
            stdev_direction_fp = statistics.stdev(all_direction_fp) 
            stdev_direction_fn = statistics.stdev(all_direction_fn) 
            stdev_d_sensitivity = statistics.stdev(all_d_sensitivity) 
            stdev_d_specificity = statistics.stdev(all_d_specificity)
            stdev_d_accuracy = statistics.stdev(all_d_accuracy)
            stdev_d_smith = statistics.stdev(all_d_smith) 

            txt_results = experiment_directory_full + "\\stdev_results.txt"
    
            with open(txt_results, 'x') as f:
                f.write("experiment: {}\n".format(subdirectory))
                f.write("stdev_connection_tp : {}\n".format(stdev_connection_tp))
                f.write("stdev_connection_tn : {}\n".format(stdev_connection_tn))
                f.write("stdev_connection_fp : {}\n".format(stdev_connection_fp))
                f.write("stdev_connection_fn : {}\n".format(stdev_connection_fn))
                f.write("stdev_c_sensitivity : {}\n".format(stdev_c_sensitivity))
                f.write("stdev_c_specificity : {}\n".format(stdev_c_specificity))
                f.write("stdev_c_accuracy : {}\n".format(stdev_c_accuracy))
                f.write("stdev_direction_tp : {}\n".format(stdev_direction_tp))
                f.write("stdev_direction_tn : {}\n".format(stdev_direction_tn))
                f.write("stdev_direction_fp : {}\n".format(stdev_direction_fp))
                f.write("stdev_direction_fn : {}\n".format(stdev_direction_fn))
                f.write("stdev_d_sensitivity : {}\n".format(stdev_d_sensitivity))
                f.write("stdev_d_specificity : {}\n".format(stdev_d_specificity))
                f.write("stdev_d_accuracy : {}\n".format(stdev_d_accuracy))
                f.write("stdev_d_smith : {}\n".format(stdev_d_smith))            
                f.close()
                
            progressbar(exp_count, NUM_EXPERIMENTS, 50, subdirectory, run)
            exp_count += 1

