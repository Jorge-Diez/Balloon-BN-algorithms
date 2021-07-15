import os
import numpy
from natsort import natsorted
from distutils.dir_util import copy_tree
from scipy import stats
import pandas as pd
import numpy as np
import statistics
import pickle



#pickle data
METRIC_SAVE_LOCATION = r"D:\My Documents\MUIA\TFM\CasualFMRI\METRICS\iMAGES_noz\Boxcar_Empirical_Resting_State\Boxcar_Empirical_Resting-State_5min\run1"


with open( METRIC_SAVE_LOCATION + "\\mean_metrics.pickle" , 'rb') as handle:
    DATA = (pickle.load(handle))