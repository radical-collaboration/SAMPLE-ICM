"""
Created on Thu Feb 15 11:10:56 2018, added to Amarel on May 14, 2018

@author: dgilford
"""

# packages needed
import numpy as np
import netCDF4 as nc
import os
import pandas as pd

## FUNCTION TO READ IN FORT.22 FILE:
# ------------------------------------------------------------------
# set up the function
def read_fort22(filepath):
    
    # get the file location
    file_loc=filepath
    
    # import the packages we need
    import numpy as np
    
    # load in the data:
    
    # get the file id
    fid = open(file_loc,'r')
    # read all of the lines into an object
    contents = fid.readlines()
    
    # intialize the numpy data array
    nvars=28
    #data=np.zeros((len(contents)-4,nvars),dtype='f')
    
    # loop over the contents, initialize
    count=0
    
    for i in range(len(contents)):
        # read in a line at a time, splitting it out into its elements
        data_split=contents[i].split()

        # figure out if we are dealing with data or header
        try:
            # this checks whether this can be turned into a number
            check=([float(ii) for ii in data_split])
            # if it can, then on the first iteration make a floating array
            if count==0:
                data_list = [[float(ii) for ii in data_split]]
                count=count+1
            # data in list form, appending created list above
            else:
                data_list.append([float(ii) for ii in data_split])
            
        except ValueError:
            if i==0:
                hline = [[str(ii) for ii in data_split]]
            else:
                hline.append([str(ii) for ii in data_split])
        
    # convert our list to a numpy array
    data_arr=np.asarray(data_list,dtype='f')
                
    # find the location of the column names
    for counter, vals in enumerate(hline):
        if vals[0]=='time':
            columnname_i=counter
            
    # create a dictionary from the headers and data
    for i in range(len(hline[columnname_i])):

        # on the first variable created the dictionary
        if i==0:
            data_out={str(hline[columnname_i][i]):data_arr[:,i]}
        # in the following variables further populate it
        else:
            data_out[str(hline[columnname_i][i])] = data_arr[:,i]
            
    # close the file id
    fid.close()            

    # return output to the above program level
    return(data_out)


## FUNCTION TO READ IN FORT.92.NC FILE:
# ------------------------------------------------------------------

# set up the function, with a variable number of calls
def read_fort92(filepath,*args):

    # set the filename and get the file's location
    file_loc=filepath
    
    # import the packages we need
    import netCDF4 as nc
    import numpy as np
    
    # load in the data:
    
    # get the file id
    fid = nc.Dataset(file_loc, mode='r')
    
    # load and store the grids as a dictionary
    grid_labels=['x0','x1','y0','y1','time']
    steps={'dx(km)':abs(fid.variables['x1'][0]-fid.variables['x1'][1]),
           'dy(km)':abs(fid.variables['y1'][0]-fid.variables['y1'][1]),
           'dt(years)':abs(fid.variables['time'][-1]-fid.variables['time'][-2])}
    grids={'x0':fid.variables['x0'][:],'x1':fid.variables['x1'][:],
           'y0':fid.variables['y0'][:],'y1':fid.variables['y1'][:],
           'time':fid.variables['time'][:]}
    
    # for each additional argument, add it to the results dictionary
    # if at least one variable is requested, loop
    if len(args)>0:
        # for all arguments (range used to make this iterable), load the data
        for i in range(len(args)):
            # if this is the first argument called, try to make the dictionary
            if i==0:
                # check if the key requested exists, if it does then store the data
                try:
                    data_out={str(args[i]):fid.variables[args[i]][:]}
                # if the key doesn't exist, print out all keys
                except KeyError:
                    print('The variable name "'+args[i]+'" does not exist, please update string...')
                    print('Available names are:')
                    # get the keys from the dictionary and loop over them plotting each
                    # note here that the list is iterable
                    for key in fid.variables.keys():
                        print(key)
                    # make sure to return so that the function doesn't continue
                    return None
            # for all others after the first, continue to populate the dictionary
            else:
                # check if the key requested exists, if it does then store the data
                try:
                    data_out[str(args[i])] = fid.variables[args[i]][:]
                # if the key doesn't exist, print out all keys
                except KeyError:
                    print('The variable name "'+args[i]+'" does not exist, please update string...')
                    print('Available names are:')
                    # get the keys from the dictionary and loop over them plotting each
                    # note here that the list is iterable
                    for key in fid.variables.keys():
                        print(key)
                    # make sure to return so that the function doesn't continue
                    return None
                
    # if no variables were requested, return "None"
    else: data_out=None
    
    # close the file id
    fid.close()
    
    # return output to the above program level
    return(grids,steps,data_out)
