library(dbnR)

data = read.csv("D:/My Documents/MUIA/TFM/CasualFMRI/NETSIM_BOLD_CENTERED/Resting_State/Resting-State_5min/run1/concatenated_bold.csv", fileEncoding="UTF-8", sep = ",", dec = ",", header = FALSE)


size <- 3
dt_train <- data[0:12750]
dt_val <- motor[12751:15000]
net <- learn_dbn_struc(dt_train, size)