# Analyze all of the data in a folder with the GES algorithm
import os
import numpy
from natsort import natsorted
from distutils.dir_util import copy_tree
from scipy import stats
import pandas as pd
import numpy as np
import statistics

NUM_RUNS= 880

#where is all the bold data
DATA_FOLDER = r"D:\My Documents\MUIA\TFM\CasualFMRI\NETSIM_Balloon_Data"

#Where BOLD data will be copied
DATA_SAVE_LOCATION = r"D:\My Documents\MUIA\TFM\CasualFMRI\NETSIM_BOLD_CENTERED"




def progressbar(acc, total, total_bar, exp, run):

    
    
    frac = acc/total
    filled_progbar = round(frac*total_bar)
    print('\r', '#'*filled_progbar + '-'*(total_bar-filled_progbar) + str(exp) + " " + str(run), '[{:>7.2%}]'.format(frac), end='')



run_count = 0

for (directory) in natsorted(os.listdir(DATA_FOLDER)):

    main_experiment_directory = os.path.join(DATA_FOLDER, directory)
    
    for (subdirectory) in natsorted(os.listdir(main_experiment_directory)):
        
        experiment_directory = os.path.join(main_experiment_directory, subdirectory)
        
        for (run) in natsorted(os.listdir(experiment_directory)):
            
            experiment_directory_run = os.path.join(experiment_directory, run)
            bold_data_directory = os.path.join(experiment_directory_run, "BOLD_DATA\FILTERED")
            
            
            direction_new_bold = DATA_SAVE_LOCATION + "\\" + directory +  "\\" + subdirectory + "\\" + run
            
            os.makedirs(direction_new_bold)
            
            all_zscored_data = np.array([])
            
            for (subject) in natsorted(os.listdir(bold_data_directory)):
                subject_data = os.path.join(bold_data_directory, subject)
                
                data = pd.read_csv(subject_data, header=None) 
                
                zscored_data_subject = stats.zscore(data, axis=0)
                
                all_zscored_data = np.vstack((all_zscored_data, zscored_data_subject)) if all_zscored_data.size else zscored_data_subject
                
                
            
            place_to_save = direction_new_bold + "\\concatenated_bold.csv"
                
            pd.DataFrame(all_zscored_data).to_csv(place_to_save, header=False, index=False)
    
            progressbar(run_count, NUM_RUNS, 50, subdirectory, run)
            run_count = run_count + 1