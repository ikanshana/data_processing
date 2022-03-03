import netCDF4 as nc
import matplotlib as mat
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import linalg as LA
import time as timer
import numpy.ma as ma

data_dir = '/workdir/indra.kanshana/2bar_project/data/'

print('POD started for combined runs')

start = timer.time()

print("\n Reading data... ")

#For 4 Hz data
fn = data_dir + 'G_H_1.25_large_run_combined.nc'

#fn = data_dir + 'combined_runs.nc'
#'/home/administrateur/Documents/lmfl_data/G_H_2.4_100Hz/combined_runs.nc'

time_filter = 1

ds = nc.Dataset(fn)

velx = ds['vel_x'][::time_filter,:,:]
vely = ds['vel_y'][::time_filter,:,:]
print("Shape of read u-component of velocity", velx.shape)
grid_x = ds['grid_x'][:,:]
grid_y = ds['grid_y'][:,:]

#time = ds['time'][::time_filter]

ds.close()

dimt, dimy, dimx = velx.shape[0], velx.shape[1], velx.shape[2]

velx_fluc = velx - np.mean(velx, axis=0, keepdims=True)
vely_fluc = vely - np.mean(vely, axis=0, keepdims=True)

velx_col = np.moveaxis(velx_fluc, 0, 2)
vely_col = np.moveaxis(vely_fluc, 0, 2)

velx_col = velx_col.reshape(dimx*dimy, dimt)
vely_col = vely_col.reshape(dimx*dimy, dimt)

U_col = np.concatenate((velx_col, vely_col), axis=0)

U_col = ma.getdata(U_col)

U_col = U_col.astype('float64')

end = timer.time()
print("\n Finished in ", end - start)


start = timer.time()
print("\n Computing Correlation matrix ...")
C = np.matmul(U_col.T, U_col)
end = timer.time()
print("\n Finished in ", end - start)


for i in range(C.shape[0]):
    for j in range(i):
        C[i,j] = C[j,i]

start = timer.time()
print("\n Solving eigenvalue problem ...")
nmodes = 10

w, v = LA.eigh(C,eigvals=(dimt - nmodes, dimt-1))
w = w[::-1]
v = v[:,::-1]

end = timer.time() 
print("\n Finished in ", end - start)


start = timer.time()
print("\n Computing POD modes ...")
U_pod =  np.matmul(U_col,v[:,:nmodes])

for i in range(nmodes):
    U_pod[:,i] = U_pod[:,i]/np.linalg.norm(U_pod[:,i])

end = timer.time()
print("\n Finished in ", end - start)


start = timer.time()
print("\n Saving POD modes ...")
fn = data_dir + 'combined_G_H_1.25_pod_modes.nc'
#'/home/administrateur/Documents/lmfl_data/G_H_2.4_100Hz/POD1/combined_' + 'k_pod_modes.nc'

try: 
    ds_w.close()  # just to be safe, make sure dataset is not already open.
except: 
    pass

ds_w = nc.Dataset(fn, 'w', format='NETCDF4')

dim_x = ds_w.createDimension('dim_x', dimx)     
dim_y = ds_w.createDimension('dim_y', dimy)    
dim_nmodes = ds_w.createDimension('nmodes', None)
dim_time = ds_w.createDimension('dim_t', dimt) 

for dim in ds_w.dimensions.items():
    print(dim)

grid_x_w = ds_w.createVariable('grid_x', np.float32, ('dim_y', 'dim_x'))
grid_y_w = ds_w.createVariable('grid_y', np.float32, ('dim_y', 'dim_x'))
#time_w = ds_w.createVariable('time',np.float32,('dim_t',))
eig_vec = ds_w.createVariable('eig_vec',np.float64,('dim_t','nmodes')) 
eig_val = ds_w.createVariable('eig_val',np.float64,('nmodes',)) 
mean_vel_x = ds_w.createVariable('mean_vel_x',np.float64,('dim_y','dim_x'))
mean_vel_y = ds_w.createVariable('mean_vel_y',np.float64,('dim_y','dim_x'))
pod_vel_x = ds_w.createVariable('pod_vel_x',np.float64,('dim_y','dim_x','nmodes'))
pod_vel_y = ds_w.createVariable('pod_vel_y',np.float64,('dim_y','dim_x','nmodes'))
coeff_t = ds_w.createVariable('coeff_t',np.float64,('nmodes','dim_t'))

grid_x_w[:,:] = grid_x
grid_y_w[:,:] = grid_y
#time_w[:] = time

eig_vec[:,:] = v
eig_val[:] = w

mean_vel_x[:,:] = np.mean(velx, axis=0)
mean_vel_y[:,:] = np.mean(vely, axis=0)

for i in range(nmodes):
    pod_vel_x[:,:,i] = U_pod[:dimx*dimy,i].reshape(dimy,dimx)
    pod_vel_y[:,:,i] = U_pod[dimx*dimy:,i].reshape(dimy,dimx)

coeff_t[:,:] = np.matmul(U_pod.T, U_col)

ds_w.close

end = timer.time()
print("\n Finished in ", end - start)

