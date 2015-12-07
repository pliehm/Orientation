# script to fit orientation of an emitter in an OLED
# INPUT: K_all.txt which contains Purcell Factor
# OUTPUT: orientation + least square error

import numpy as np
import time
import matplotlib.pyplot as plt
import os as os

###################
# PARAMETER INPUT #
###################

# radiative efficiency in vacuum
eta = 0.8

# give range for the wavelength used in the simulation and the measurement
wave_min = 530
wave_max = 750
wave_step = 1

# give range for the angles used in simulation and experiment
angle_min = 0
angle_max = 89
angle_step = 1

# sepcify folder where the experimental data is
exp_data_folder = "Samples"

# specify folder with the simulation for the different cavity thickness
cavities = "cavities"


#########
# START #
#########

#########################
# Define some functions #
#########################


# define a function which calculates the Purcell Factor weigthed by the anisotropy factor a
def F_func(a):
    F = (1-a)*(F_tmh + F_teh)+a*F_tmv
    return F



# define function for eta with arrays
def eta_func(a,eta): # "a" is the anisotropy factor, "et" is the non weighted radiative efficiency, e.g. measure independently 
	# calculate effective radiative efficiency (Eq.20 in M. Furnos paper )
    eta_eff=eta*F_func(a)/(1-eta+F_func(a)*eta)
    
    # do we need to normalise here?
    eta_eff = eta_eff/eta_eff.max()    
    return eta_eff

# function to calculate the normalised spectral radiant intensity 
def I_sim(a):
    # calculate the radiative efficiency
    eta_eff=eta_func(a,eta)

    # calculate spectral radiant intensity (Eq. A19 in M. Furnos paper )
    I=((3/2.0)*(1-a)*I_sim_TMh.T+3*a*I_sim_TMv.T)*eta_eff/F_func(a)/wave_sim
    I=I.T
    
    # with lists
#    for i in range(len(I_sim_TMh)):
#        I.append(((3/2.0)*(1-a)*float(I_sim_TMh[i])+3*a*float(I_sim_TMv[i]))*eta_wave[i]/F_func(a,i)/wave_sim[i])
    
    # normalise the intensity to compare with experiment
    I=I/I.max()
    return I
    
# do all the calculation in one function to make code faster

# calculate how many angles 
angles = (angle_max-angle_min)/angle_step + 1

#  create array with wavelengths
wave_sim = np.arange(wave_min,wave_max+1,wave_step, dtype=np.float)

waves = len(wave_sim)
    
# make folderlist of folder "samples"
samples = os.listdir(exp_data_folder)

# sort the list
samples.sort()

# list to store results of samples
sample_results = []
# for all folders in the folder "exp_data_folder" do:
for sample in samples[0:3]:
    
    # make list to store: cavity thickness, best fitted a, error
    result = []
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
    
    # only take the intensity of file and reshape it to a 2D array (matrix)
    # get size of the array from parameters
    I_exp = np.reshape(p[:,2],(angles+1,waves)).T
    
    # remove last angle
    I_exp = I_exp[:,:-1]
    I_exp = I_exp/I_exp.max()
    plt.imshow(I_exp, aspect=angles/float(waves))
    plt.show()
    

    for cavity in cavities_list:
         #show which thickness is fitted
        #print 'thickness:', cavity, 'nm'
        
        # load simulated spectral radiant intensity
        p = np.loadtxt(cavities +'/'+ cavity + '/SpectralRadiantIntensity_TMh.txt')
        I_sim_TMh = np.reshape(p[:,2],(waves,angles))
        #plt.imshow(I_sim_TMh, aspect=0.084)
        #plt.show()
        
        p = np.loadtxt(cavities +'/'+ cavity + '/SpectralRadiantIntensity_TMv.txt')
        I_sim_TMv = np.reshape(p[:,2],(waves,angles))
        #plt.imshow(I_sim_TMv,aspect=0.084)
        #plt.show()
        
        # load the purcell factor        
        p = np.loadtxt(cavities +'/'+ cavity + '/K_all.txt')
        F_tmh = p[:,0]
        F_tmv = p[:,1]
        F_teh = p[:,2]
        

        
        # fit anisotropy factor a
        # create list for the error
        error_list = []
        a_list = []
        for a in np.arange(0,1,0.01):

            # calculate difference between simulation and experiment for a given a
            error = ((I_exp-I_sim(a))**2).sum()
            # append a and the difference to a list
            error_list.append(error)
            a_list.append(a)
        
        # get the a with the smallest error
        a_fitted = a_list[error_list.index(min(error_list))]
        
        # store results
        result.append([int(cavity),a_fitted, min(error_list)])
               
        
        #print a_fitted, min(error_list)
    
    result = np.asarray(result)

    
    # take best result for cavity
    result_best = result[np.where(result[:,2]==result[:,2].min())]
    
    print "Best result:\n", result_best
    
    # plot a versus the cavity thickness    
    plt.figure(1)
    plt.plot(result[:,0], result[:,1], 'o')
    plt.xlabel('ETL thickness in nm')
    plt.ylabel('anisotropy factor')
    plt.grid()
    plt.savefig(str(sample)[:5] + '_results/'+'a')
    #plt.show()
    plt.clf()    
        
    # plot the error versus the cavits thickness
    plt.figure(2)
    plt.plot(result[:,0], result[:,2], 'o')
    plt.xlabel('ETL thickness in nm')
    plt.ylabel('leastsquare')
    plt.grid()
    plt.savefig(str(sample)[:5]+'_results/'+'leastsquare')
    #plt.show()
    plt.clf()    
       
       
    print "Time needed: ", time.clock()-start," seconds"

        
        
    
    