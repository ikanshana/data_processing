import netCDF4 as nc
import matplotlib as mat
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy.ma as ma
import glob
import sys
import os
import argparse
from pathlib import Path


sys.path.append(os.path.abspath(
    "/home/administrateur/Documents/lmfl_data/Programmes/python_codes/langevin-regression"))

import utils


def main():

    parser = argparse.ArgumentParser(
        description="Compute time series, drift, diffusion and pdf from .nc files.")
    parser.add_argument("datadir", type=str, help="path to the datafiles.")
    parser.add_argument("N_bins", type=int, help="Number of bins for pdf.")
    parser.add_argument("stride", type=int, help="How many snapshots to skip.")
    parser.add_argument("f_sample", type=int, help="Sampling freauency of PIV dataset.")
    parser.add_argument("y_cut_plus", type=float, help="Positive cut-off for state.")
    parser.add_argument("y_cut_minus", type=float, help="Negative cut-off for state.")

    args = vars(parser.parse_args())

    datadir = args["datadir"]
    N_bins = args["N_bins"]
    stride = args["stride"]
    f_sample = args["f_sample"]
    y_cut_plus = args["y_cut_plus"]
    y_cut_minus = args["y_cut_minus"]

    Y_m, time, total_transitions, all_transition_instants = Combined_Ym(datadir, y_cut_plus, y_cut_minus)
    Plot_Ym(datadir, time, Y_m, all_transition_instants, y_cut_plus, y_cut_minus)
    Plot_DDP(datadir, N_bins, Y_m, stride, f_sample)


def Average_Y(u,v,x,y):
    #u,v,x,y are two-dimensional fields
    K = u**2 + v**2
    Y_m = np.sum(y*K, axis=(1, 2))/np.sum(K, axis=(1, 2))
    return Y_m


def Combined_Ym(datadir, y_cut_plus, y_cut_minus):

    all_files = sorted(glob.glob(datadir + "*.nc"))

    i = 0
    total_transitions = 0
    all_transition_instants = []


    for file in all_files:
        ds = nc.Dataset(file)

        velx = ds['vel_x']
        vely = ds['vel_z']

        grid_x = ma.getdata(ds['grid_x'][:,:])
        grid_y = ma.getdata(ds['grid_z'][:,:])

        time1 = ma.getdata(ds['time'][:])
        velx = ma.getdata(velx[:,0:,0:])
        vely = ma.getdata(vely[:,0:,0:])

        if i == 0:
            N_t,N_files = len(velx), len(all_files)
            Y_m = np.zeros(N_t*N_files)
            time = np.zeros(N_t*N_files)


        Y_m1 = Average_Y(velx, vely, grid_x, grid_y)
        transition_instants, state_transition_count = switch_counts(Y_m1, time1, y_cut_plus, y_cut_minus)

        if state_transition_count != 0 :
            print("Switch for run ",i+1,  " at time instants : ", transition_instants)
            total_transitions = total_transitions + state_transition_count

            for transition_instant in transition_instants :
                all_transition_instants = all_transition_instants + [time[i*N_t - 1] + transition_instant]

        Y_m[i*N_t:(i+1)*N_t] = Y_m1
        time[i*N_t:(i+1)*N_t] = time[i*N_t - 1] + time1

        i = i + 1

    print("Total count of switches: ", total_transitions)

    return Y_m, time, total_transitions, all_transition_instants

def switch_counts(Y_m, time, y_cut_plus, y_cut_minus, series_count = 1, sample_N = 0):

    """
    If series_count is specified more than 1, it needs sample size of the series
    , sample_N to discount switches coming from jump between two consecutive series
    """

    #Initial states and counters
    #+ is left, - is right
    if Y_m[0] > y_cut_plus:
        curr_state = 'left'
    elif Y_m[0] < y_cut_minus:
        curr_state = 'right'
    else :
        curr_state = 'center'

    state_transition_count = 0
    switch = 0
    prev_state = curr_state
    prev_prev_state = curr_state
    counter = -1
    transition_instants = []


    for t, y in zip(time, Y_m):

        counter = counter + 1
        if y > y_cut_plus:
            curr_state = 'left'
        elif y < y_cut_minus:
            curr_state = 'right'
        else :
            curr_state = 'center'

        if curr_state != prev_state:
            switch = 1
        else :
            switch = 0

        if switch == 1:
            if curr_state != 'center' and curr_state != prev_prev_state and prev_prev_state != 'center':

            #If multiple series than check if the switch is due to series change
                if series_count == 1:
                    state_transition_count = state_transition_count + 1
                    transition_instants.append(t)

                else :
                    if (t/(time[1] - time[0])) != sample_N:
                        state_transition_count = state_transition_count + 1
                        transition_instants.append(t)
                    else:
                        prev_state = curr_state
                        prev_prev_state = curr_state

            prev_prev_state = prev_state

        prev_state = curr_state


    return transition_instants, state_transition_count


def Plot_Ym(datadir, time, Y_m, all_transition_instants, y_cut_plus, y_cut_minus):

    ##Time series of combined runs
    fg_Y = plt.figure(figsize=(14, 6))
    font = {'family': 'normal', 'weight': 'normal', 'size': 32}
    plt.rc('font', **font)

    ax1 = fg_Y.add_subplot(111)
    ax1.plot(time, Y_m)
    ax1.set_xlabel('time (s)')
    ax1.set_ylabel('$Y_m$')
    plt.axhline(y= y_cut_plus,ls='--',c= 'g',label = 'state cut-off')
    plt.axhline(y= y_cut_minus,ls='--',c= 'g')


    for transition_instant in all_transition_instants[:-1]:
        plt.axvline(x=transition_instant,ls='--',c= 'r')

    #Use last instant for legend
    plt.axvline(x=all_transition_instants[-1],ls='--',c= 'r', label = 'switch instant')

    plt.legend(loc='lower right', bbox_to_anchor=(1, 1),ncol=3)
    plt.tight_layout()
    plt.savefig(datadir + "Ym_series.jpg", dpi=150)
    #plt.show()

#Plot drift, diffusion and pdf distribution
def Plot_DDP(datadir, N_bins, Y_m, stride, f_sample):

    ### Compute K-M coeffs + Edges ###
    dt = 1/f_sample

    bins = np.array([N_bins])

    ## The size is only valid for 1-D case
    KMc_exp = np.zeros((2,N_bins))


    Edges_X = np.linspace(min(Y_m), max(Y_m), N_bins + 1)
    X_values = (np.diff(Edges_X, axis=0) + 2*Edges_X[:-1])/2

    ### Compute histogram ###

    p_hist = np.histogramdd(Y_m, bins=bins, density=True)[0]

    ### Build KM coeffs ###

    f_KM, a_KM, f_err, a_err = utils.KM_avg(Y_m, Edges_X, stride=stride, dt=dt)

    ### Fill with NaN for values not present in the interval ###
    filtre = p_hist == 0
    f_KM[filtre], a_KM[filtre], f_err[filtre], a_err[filtre] = np.NaN, np.NaN, np.NaN, np.NaN


    KMc_exp[0, :] = f_KM
    KMc_exp[1, :] = a_KM


    ##Plot drift, diffusion and PDF plots ###
    fg_dd_fit = plt.figure(figsize=(16, 10))
    ax1 = fg_dd_fit.add_subplot(221)
    ax2 = fg_dd_fit.add_subplot(222)
    ax3 = fg_dd_fit.add_subplot(223)

    font = {'family': 'normal', 'weight': 'normal', 'size': 32}

    plt.rc('font', **font)

    ax1.plot(X_values, f_KM, 'ko')
    ax1.set_xlabel('Y')
    ax1.set_ylabel('f(Y)')
    ax1.ticklabel_format(axis='y', style='sci', scilimits=(-1,1))

    ax2.plot(X_values, a_KM, 'ko')
    ax2.set_xlabel('Y')
    ax2.set_ylabel('a(Y)')
    ax2.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))

    ax3.plot(X_values, p_hist, 'ko')
    ax3.set_xlabel('Y')
    ax3.set_ylabel('p(Y)')

    plt.tight_layout()
    fg_dd_fit.subplots_adjust(wspace=0.2, hspace=0.4, top=0.85)

    fg_dd_fit.suptitle("Drift, diffusion and PDF plots")

    plt.savefig(datadir + "DDP.jpg", dpi=150)

    #plt.show()


if __name__ == "__main__":
    main()
