# Analyze all of the data in a folder with the GES algorithm
import os
import numpy
from natsort import natsorted

NUM_EXPERIMENTS = 52800

#where is all the bold data
DATA_FOLDER = r"/home/jorge/Desktop/CasualFMRI/NETSIM_Balloon_Data"

#Where BOLD data will be copied
DATA_SAVE_LOCATION = r"/home/jorge/Desktop/CasualFMRI/NETSIM_ALGORITHM_RESULTS/FGES"

#Where our java is located, and other parameters
ccd_loc =  r'"/home/jorge/Desktop/CasualFMRI/CasualCCD/causal-cmd-1.2.2/causal-cmd-1.2.2-jar-with-dependencies.jar"'
algorithm = "--algorithm fges --penaltyDiscount 1 --semBicStructurePrior 1"
before_data = "--data-type continuous --dataset"

csv_options_out = "--delimiter comma --no-header --out"
score_verbose = "--score sem-bic-score --verbose"



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
            
            direction_new_bold = DATA_SAVE_LOCATION + "/" + directory +  "/" + subdirectory + "/" + run 
            
            if (not (os.path.isdir(direction_new_bold))):
                os.makedirs(direction_new_bold)
            
            
            
            for (subject) in natsorted(os.listdir(bold_data_directory)):
                subject_data = os.path.join(bold_data_directory, subject)
                exp_count = exp_count + 1
                
                sub_out = subject.rsplit("bold_filtered")[0]
                
                command = "java -jar " + ccd_loc + " " + algorithm + " " + before_data + " " + "\"" + subject_data  + "\"" + \
                     " " + csv_options_out + " " + "\"" + direction_new_bold + "\"" + " --prefix " + sub_out + " "   +  " " + score_verbose + \
                     "> /dev/null"

                
                #FALTA --PREFIX DESPUES DE EL OUT, ESTE SERA EL NOMBRE
                
                
                
                os.system(command)
                progressbar(exp_count, NUM_EXPERIMENTS, 50, subdirectory, run)