import netCDF4 as nc
import matplotlib as mat
import numpy as np
import matplotlib.pyplot as plt
import os
from numpy import linalg as LA
import time as timer

data_dir = '/home/administrateur/Documents/lmfl_data/G_H_2.4/'

i = 0
for run in range(1,11):

    print('Reading run ' + str(run))

    if run < 10:
        fn = data_dir + 'G_H_2.4_large_run0' + str(run) + '.nc'
    else:
        fn = data_dir + 'G_H_2.4_large_run' + str(run) + '.nc'

    ds = nc.Dataset(fn)

    velx = ds['vel_x'][:,:,:]
    vely = ds['vel_y'][:,:,:]

    grid_x = ds['grid_x'][:,:]
    grid_y = ds['grid_y'][:,:]
    flag = ds['flag'][:,:]
    dim_t =  ds.dimensions['dim_t'].size

    #time = ds['time'][:]
    ds.close()


    #taking subset of the image for the POD
    velx = velx[:, np.any(abs(grid_y)< 0.045, axis = 1), :][:,:, np.any(abs(grid_x)< 0.25, axis = 0)]
    vely = vely[:, np.any(abs(grid_y)< 0.045, axis = 1), :][:,:, np.any(abs(grid_x)< 0.25, axis = 0)]
    flag = flag[:, np.any(abs(grid_y)< 0.045, axis = 1), :][:,:, np.any(abs(grid_x)< 0.25, axis = 0)]


    grid_x1 = grid_x[np.any(abs(grid_y)< 0.045, axis = 1), :][:, np.any(abs(grid_x)< 0.25, axis = 0)]
    grid_y = grid_y[np.any(abs(grid_y)< 0.045, axis = 1), :][:, np.any(abs(grid_x)< 0.25, axis = 0)]

    grid_x = grid_x1

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
        grid_y_w = ds_w.createVariable('grid_y', np.float32, ('dim_y', 'dim_x'))
        grid_x_w[:,:] = grid_x
        grid_y_w[:,:] = grid_y

        vel_x_w = ds_w.createVariable('vel_x', np.float64,('dim_t', 'dim_y','dim_x'))
        vel_y_w = ds_w.createVariable('vel_y', np.float64,('dim_t', 'dim_y','dim_x'))
        flag_w = ds_w.createVariable('flag', np.float64,('dim_t', 'dim_y','dim_x'))
        #time_w = ds_w.createVariable('time',np.float32,('dim_t',))

    #####  Add for the ith iteration to combined data
    vel_x_w[i*dim_t:(i+1)*dim_t,:,:] = velx
    vel_y_w[i*dim_t:(i+1)*dim_t,:,:] = vely
    flag_w[i*dim_t:(i+1)*dim_t,:,:] = flag
    #time_w[i*dim_t:(i+1)*dim_t] = time
    i = i + 1
    print('Added data for run ' + str(run))

ds_w.close()
