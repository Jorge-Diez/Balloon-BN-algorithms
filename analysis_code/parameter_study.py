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

metric_data_origin = r"D:\My Documents\MUIA\TFM\CasualFMRI\METRICS_NOCYCLE_PARAM"

sensitivity_c = []
specificity_c = []

sensitivity_d = []
specificity_d = []


algorithms = ["DIRECTLINGAM", "FGES", "GOBLIN", "IMAGES", "LINGAM", "MMHC", "TABU" ]
exp_names = []


structures = ["Resting-State", "Boxcar_Resting-State", "Empirical-Resting-State", "Boxcar_Empirical_Resting-State"]

structures_simple = ["5nodeR", "5nodeB", "8nodeR", "8nodeB"]

string_experiments = ["60sub-1200ms-10min", "5min", "30min", "60min", "10sub", "30sub", "100sub", "tr250ms", "tr750ms", "tr3000ms", "excellentconditions"]


initial_point = r"D:\My Documents\MUIA\TFM\CasualFMRI\METRICS_NOCYCLE_PARAM\DIRECTLINGAM"





for i,structure in enumerate(structures):
    
    for param_exp in string_experiments:
        
        exp_full = structure + "_" + param_exp
        
        exp_to_join = os.path.join(structure, exp_full)
        
        line_sensitivity_c = []
        line_specificity_c = []
        
        line_sensitivity_d = []
        line_specificity_d = []
        
        exp_names.append(structures_simple[i] + " " + param_exp)
        
        
        for algorithm in algorithms:
            algorithm_direction = os.path.join(metric_data_origin, algorithm)
            
            
            
            full_path_param = os.path.join(algorithm_direction, exp_to_join)
            
            
            data_location_txt = full_path_param + "\\results.txt"
            
            file = open(data_location_txt)
            metrics = file.readlines()
            
            sensitivity_line_c = metrics[5]
            specificity_line_c = metrics[6]
            
            sensitivity_line_d = metrics[12]
            specificity_line_d = metrics[13]
            
            
            
            line_sensitivity_c.append(float((sensitivity_line_c.rsplit(": ")[1])[:-2]))
            line_specificity_c.append(float((specificity_line_c.rsplit(": ")[1])[:-2]))
            
            line_sensitivity_d.append(float((sensitivity_line_d.rsplit(": ")[1])[:-2]))
            line_specificity_d.append(float((specificity_line_d.rsplit(": ")[1])[:-2]))
                    
        sensitivity_c.append(line_sensitivity_c) 
        specificity_c.append(line_specificity_c)  
        
        sensitivity_d.append(line_sensitivity_d) 
        specificity_d.append(line_specificity_d) 
        
        

sensitivity_c = np.array(sensitivity_c)
specificity_c = np.array(specificity_c)

sensitivity_d = np.array(sensitivity_d)
specificity_d = np.array(specificity_d)








plt.title("specificity-d", fontsize = 32)







fig = plt.figure(figsize=(7,44))
ax = fig.add_subplot(111)

cax = ax.matshow(specificity_d, cmap = 'plasma', vmin=0,vmax=1)

ax.set_xticklabels(['']+algorithms)
ax.set_yticks(np.arange(0, 44, 1))
ax.set_yticklabels(exp_names)


plt.xticks(rotation=90)

ax.plot(  [-0.5,6.5] , [10.5,10.5], 'k'   )
ax.plot(  [-0.5,6.5] , [21.5,21.5], 'k'   )
ax.plot(  [-0.5,6.5] , [32.5,32.5], 'k'   )


fig.colorbar(cax)



plt.show()






