# Analyze all of the data in a folder with the GES algorithm
import os
import numpy
from natsort import natsorted

NUM_EXPERIMENTS = 880

#where is all the bold data
DATA_FOLDER = r"D:\My Documents\MUIA\TFM\CasualFMRI\NETSIM_BOLD_CENTERED"

#Where BOLD data will be copied
DATA_SAVE_LOCATION = r"D:\My Documents\MUIA\TFM\CasualFMRI\NETSIM_ALGORITHM_RESULTS\iMAGES"

#Where our java is located, and other parameters
ccd_loc =  r'"D:\My Documents\MUIA\TFM\CasualFMRI\causal-cmd-1.3.0/causal-cmd-1.3.0-jar-with-dependencies.jar"'
algorithm = "--algorithm imgs_cont --faithfulnessAssumed 1 --semBicStructurePrior 1"
before_data = "--data-type continuous --dataset"

csv_options_out = "--delimiter comma --no-header --out"
score_verbose = "--verbose"



def progressbar(acc, total, total_bar, exp, run):

    
    
    frac = acc/total
    filled_progbar = round(frac*total_bar)
    print('\r', '#'*filled_progbar + '-'*(total_bar-filled_progbar) + str(exp) + "  " +  str(run) + "|| Subject count: " + str(acc), '[{:>7.2%}]'.format(frac), end='')


exp_count = 0

for (directory) in natsorted(os.listdir(DATA_FOLDER)):

    main_experiment_directory = os.path.join(DATA_FOLDER, directory)
    
    for (subdirectory) in natsorted(os.listdir(main_experiment_directory)):
        
        experiment_directory = os.path.join(main_experiment_directory, subdirectory)
        
        for (run) in natsorted(os.listdir(experiment_directory)):
            
            experiment_directory_run = os.path.join(experiment_directory, run)
            
            bold_data_directory = os.path.join(experiment_directory_run, "BOLD_DATA/FILTERED")
            
            direction_new_bold = DATA_SAVE_LOCATION + "\\" + directory +  "\\" + subdirectory + "\\" + run 
            
            if (not (os.path.isdir(direction_new_bold))):
                os.makedirs(direction_new_bold)
                
                concatenated_data = os.path.join(experiment_directory_run, "concatenated_bold.csv")
                
                sub_out = "concatenated_results"
                command = "java -jar " + ccd_loc + " " + algorithm + " " + before_data + " " + "\"" + concatenated_data  + "\"" + \
                     " " + csv_options_out + " " + "\"" + direction_new_bold + "\"" + " --prefix " + sub_out + " "   +  " " + score_verbose + \
                         " /f > nul"
            
            
            

                
                
                exp_count = exp_count + 1
                os.system(command)
                progressbar(exp_count, NUM_EXPERIMENTS, 50, subdirectory, run)