import netCDF4 as nc
import matplotlib as mat
import numpy as np
import matplotlib.pyplot as plt
import os
from numpy import linalg as LA
import time as timer

data_dir = '/workdir/indra.kanshana/2bar_project/data/'

i = 0
for run in [3,16,17,18,19,20,21,26,35,52,103,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125]:
    
    print('Reading run ' + str(run))
    
    if run < 10:
        fn = data_dir + 'V_5mps_G_H_2.4_100Hz_run_0' + str(run) + '.nc'
    else:
        fn = data_dir + 'V_5mps_G_H_2.4_100Hz_run_' + str(run) + '.nc'

    ds = nc.Dataset(fn)
    
    velx = ds['vel_x'][:,:,:]
    velz = ds['vel_z'][:,:,:]

    grid_x = ds['grid_x'][:,:]
    grid_z = ds['grid_z'][:,:]
    flag = ds['flag'][:,:]

    time = ds['time'][:]
    ds.close()


    dimt, dimy, dimx = velx.shape[0], velx.shape[1], velx.shape[2]
    
    #####  Finished reading #####
 

    #####  Create an nc file to write for first run and add to it later#####
    if i==0:
        fn = data_dir + 'combined_runs.nc'

        try: 
            ds_w.close()  # just to be safe, make sure dataset is not already open.
        except: 
            pass
        ds_w = nc.Dataset(fn, 'w', format='NETCDF4')

        dim_time = ds_w.createDimension('dim_t', None) 
        dim_x = ds_w.createDimension('dim_x', dimx)     
        dim_y = ds_w.createDimension('dim_y', dimy)    

        grid_x_w = ds_w.createVariable('grid_x', np.float32, ('dim_y', 'dim_x'))
        grid_z_w = ds_w.createVariable('grid_z', np.float32, ('dim_y', 'dim_x'))
        grid_x_w[:,:] = grid_x
        grid_z_w[:,:] = grid_z

        vel_x_w = ds_w.createVariable('vel_x', np.float64,('dim_t', 'dim_y','dim_x'))
        vel_z_w = ds_w.createVariable('vel_z', np.float64,('dim_t', 'dim_y','dim_x'))
        flag_w = ds_w.createVariable('flag', np.float64,('dim_t', 'dim_y','dim_x'))
        time_w = ds_w.createVariable('time',np.float32,('dim_t',))

    #####  Add for the ith iteration to combined data    
    vel_x_w[i*len(time):(i+1)*len(time),:,:] = velx
    vel_z_w[i*len(time):(i+1)*len(time),:,:] = velz
    flag_w[i*len(time):(i+1)*len(time),:,:] = flag
    time_w[i*len(time):(i+1)*len(time)] = time
    i = i + 1    
    print('Added data for run ' + str(run))

ds_w.close()


