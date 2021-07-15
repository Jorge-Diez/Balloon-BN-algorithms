import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
import os
from natsort import natsorted
from distutils.dir_util import copy_tree
from scipy import stats
import pandas as pd
import numpy as np
import statistics

#sensitivity = recall
#

metric_data_origin = r"D:\My Documents\MUIA\TFM\CasualFMRI\METRICS_NOCYCLE"


title = "ILP ALGORITHM ON RESTING STATE"

labels_ALGORITHMS = ["DIRECTLINGAM", "FGES", "GOBNILP", "iMAGES", "ICA-LINGAM", "MMHC", "TABU" ]


sensitivity_c = []
specificity_c = []


sensitivity_d = []
specificity_d = []


experiment_name = "Boxcar_Resting_State"
subtype = "Boxcar_Resting-State_excellentconditions"
string_experiments = ["60sub-1200ms-10min", "5min", "30min", "60min", "10sub", "30sub", "100sub", "tr250ms", "tr750ms", "tr3000ms", "excellentconditions"]



for (algorithm) in natsorted(os.listdir(metric_data_origin)):

    main_algorithm_directory = os.path.join(metric_data_origin, algorithm)
    
    for (subdirectory) in natsorted(os.listdir(main_algorithm_directory)):
        
        if (subdirectory == experiment_name):
            
            experiment_directory = os.path.join(main_algorithm_directory, subdirectory)
        
            for (parameter_experiment) in natsorted(os.listdir(experiment_directory)):
            
                if (parameter_experiment == subtype):
                
                    parameter_experiment_directory = os.path.join(experiment_directory, parameter_experiment)
                    
                    
                    data_location_txt = parameter_experiment_directory + "//results.txt"
    
                    file = open(data_location_txt)
                    metrics = file.readlines()
                    
                    sensitivity_line_c = metrics[5]
                    specificity_line_c = metrics[6]
                    
                    sensitivity_line_d = metrics[12]
                    specificity_line_d = metrics[13]
                    
                    
                    
                    sensitivity_c.append(float((sensitivity_line_c.rsplit(": ")[1])[:-2]))
                    specificity_c.append(float((specificity_line_c.rsplit(": ")[1])[:-2]))
                    
                    sensitivity_d.append(float((sensitivity_line_d.rsplit(": ")[1])[:-2]))
                    specificity_d.append(float((specificity_line_d.rsplit(": ")[1])[:-2]))
                    
                    
         
plt.rcParams['xtick.labelsize'] = 14 
plt.rcParams['ytick.labelsize'] = 12            
            
x = np.arange(len(labels_ALGORITHMS))  # the label locations
width = 0.35  # the width of the bars                  
                    
                    
fig, ax = plt.subplots(figsize=(15,7))
rects1 = ax.bar(x - width/2, sensitivity_c, width, label='sensitivity')
rects2 = ax.bar(x + width/2, specificity_c, width, label='specificity')


ax.set_ylabel('Scores', fontsize = 14)
ax.set_title("Connections in excellent conditions, 5-node Boxcar", fontsize = 18)
ax.set_xticks(x)
ax.set_xticklabels(labels_ALGORITHMS)
ax.legend()



#ax.bar_label(rects1, padding=3)
#ax.bar_label(rects2, padding=3)
fig.tight_layout()
plt.savefig("smith_boxcar_excellent_normal_connection" , bbox_inches='tight')

plt.show()


fig2, ax2 = plt.subplots(figsize=(15,7))
rects1 = ax2.bar(x - width/2, sensitivity_d, width, label='sensitivity')
rects2 = ax2.bar(x + width/2, specificity_d, width, label='specificity')


ax2.set_ylabel('Scores', fontsize = 14)
ax2.set_title("Orientations in excellent conditions, 5-node Boxcar", fontsize = 18)
ax2.set_xticks(x)
ax2.set_xticklabels(labels_ALGORITHMS)
ax2.legend()

#ax2.bar_label(rects1, padding=3)
#ax2.bar_label(rects2, padding=3)
fig2.tight_layout()
plt.savefig("smith_boxcar_excellent_normal_orientation" , bbox_inches='tight')


plt.show()
            
            
            
            
            
            






