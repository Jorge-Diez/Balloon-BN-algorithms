import matplotlib.pyplot as plt
import numpy as np
import os
import pickle

#sensitivity = recall
#

metric_data_origin = r"D:\My Documents\MUIA\TFM\CasualFMRI\METRICS\MMHC\Resting_State"


title = "ILP ALGORITHM ON RESTING STATE"

labels = ["normal", "5min", "30min", "60min", "10sub", "30sub", "100sub", "tr250", "tr750", "tr3000", "sota" ]


sensitivity_c = []
specificity_c = []


sensitivity_d = []
specificity_d = []


experiment_name = "Resting-state_"
string_experiments = ["60sub-1200ms-10min", "5min", "30min", "60min", "10sub", "30sub", "100sub", "tr250ms", "tr750ms", "tr3000ms", "excellentconditions"]


for string_subexperiment in string_experiments:
    full_name = experiment_name + string_subexperiment
    
    data_location_txt = metric_data_origin + "\\" + full_name + "\\results.txt"
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
    
    
x = np.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars    



fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, sensitivity_c, width, label='sensitivity_c')
rects2 = ax.bar(x + width/2, specificity_c, width, label='specificity_c')


ax.set_ylabel('Scores in range 0-1')
ax.set_title("FGES connection data")
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()


ax.bar_label(rects1, padding=3)
ax.bar_label(rects2, padding=3)
fig.tight_layout()




fig2, ax2 = plt.subplots()
rects1 = ax2.bar(x - width/2, sensitivity_d, width, label='sensitivity_c')
rects2 = ax2.bar(x + width/2, specificity_d, width, label='specificity_c')


ax2.set_ylabel('Scores in range 0-1')
ax2.set_title("FGES direction data")
ax2.set_xticks(x)
ax2.set_xticklabels(labels)
ax2.legend()


ax2.bar_label(rects1, padding=3)
ax2.bar_label(rects2, padding=3)
fig2.tight_layout()

plt.show()