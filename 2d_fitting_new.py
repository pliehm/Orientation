# script to fit orientation of an emitter in an OLED
# INPUT: K_all.txt which contains Purcell Factor
# OUTPUT: orientation + least square error

import numpy as np
import time
import matplotlib.pyplot as plt
import os as os

eta = 0.8
min_wave = 475
max_wave = 700

min_angle = 0
max_angle = 85
step_angle = 5

# sepcify folder where the experimental data is
exp_data_folder = "samples"

# specify folder with the simulation for the different cavity thickness
cavities = "cavities"


#########
# START #
#########

# make folderlist of folder "samples"
samples = os.listdir(exp_data_folder)

# sort the list
samples.sort()

# for all folders in the folder "exp_data_folder" do:
for sample in samples:
    
    # print to the user part of the folder name
    print str(sample)[:5] + '\n'
    
    # make a list of all folders of the simulated cavity thicknesses
    cavities_list = os.listdir(cavities)
    # sort the list
    cavities_list.sort()
    
    # make a new folder for the results
    os.mkdir(str(sample)[:5]+'_results')
    
    #start a timer to measure how long the calculation takes, not important for the actual calculation
    start = time.clock() 

    ##########################    
    # read experimental data #
    ##########################

    print 'read experimental data'
    
    # use numpy to load the file
    p = np.loadtxt(exp_data_folder + '/' + sample,skiprows=1)
    
    # only take the intensity of file and reshape it to a 2D array
    reshaped = np.reshape(p[:,2],(226,19)).T
    
    
    
    