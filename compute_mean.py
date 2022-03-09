import netCDF4 as nc
import matplotlib as mat
import numpy as np
import matplotlib.pyplot as plt
import os
from numpy import linalg as LA
import time as timer

data_dir = '/home/administrateur/Documents/lmfl_data/G_H_2.4_100Hz/'

i = 0
run_list = [3,16,17,18,19,20,21,26,35,52,103,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125]
N_run = len(run_list)

for run in run_list:

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
    ds.close()


    dimt, dimz, dimx = velx.shape[0], velx.shape[1], velx.shape[2]

    #####  Finished reading #####


    #####  Create an nc file to write and initialize sum for average in the first run  #####
    if i==0:


        fn = data_dir + 'mean.nc'

        try:
            ds_w.close()  # just to be safe, make sure dataset is not already open.
        except:
            pass
        ds_w = nc.Dataset(fn, 'w', format='NETCDF4')

        dim_time = ds_w.createDimension('dim_t', None)
        dim_x = ds_w.createDimension('dim_x', dimx)
        dim_z = ds_w.createDimension('dim_z', dimz)

        grid_x_w = ds_w.createVariable('grid_x', np.float32, ('dim_z', 'dim_x'))
        grid_z_w = ds_w.createVariable('grid_z', np.float32, ('dim_z', 'dim_x'))
        grid_x_w[:,:] = grid_x
        grid_z_w[:,:] = grid_z

        vel_x_mean = ds_w.createVariable('vel_x', np.float64,('dim_t', 'dim_z','dim_x'))
        vel_z_mean = ds_w.createVariable('vel_z', np.float64,('dim_t', 'dim_z','dim_x'))

        vel_x_mean[0,:,:] = np.mean(velx, axis=0, keepdims=True)
        vel_z_mean[0,:,:] = np.mean(velz, axis=0, keepdims=True)

    else:
        vel_x_mean[0,:,:] = vel_x_mean[0,:,:] + np.mean(velx, axis=0, keepdims=True)
        vel_z_mean[0,:,:] = vel_z_mean[0,:,:] + np.mean(velz, axis=0, keepdims=True)

    #####  Add for the ith iteration to combined data
    i = i + 1
    print('Added data for run ' + str(run))

vel_x_mean[0,:,:] = vel_x_mean[0,:,:]/N_run
vel_z_mean[0,:,:] = vel_z_mean[0,:,:]/N_run

ds_w.close()
