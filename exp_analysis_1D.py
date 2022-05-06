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
    parser.add_argument("f_sample", type=int,
                        help="Sampling freauency of PIV dataset.")
    parser.add_argument("y_cut_plus", type=float,
                        help="Positive cut-off for state.")
    parser.add_argument("y_cut_minus", type=float,
                        help="Negative cut-off for state.")

    args = vars(parser.parse_args())

    datadir = args["datadir"]
    N_bins = args["N_bins"]
    stride = args["stride"]
    f_sample = args["f_sample"]
    y_cut_plus = args["y_cut_plus"]
    y_cut_minus = args["y_cut_minus"]

    Y_m, time, total_transitions, all_transition_instants = Combined_Ym(
        datadir, y_cut_plus, y_cut_minus)
    Plot_Ym(datadir, time, Y_m, all_transition_instants, y_cut_plus, y_cut_minus)
    Plot_DDP(datadir, N_bins, Y_m, stride, f_sample)


def Average_Y(u, v, x, y):
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

        grid_x = ma.getdata(ds['grid_x'][:, :])
        grid_y = ma.getdata(ds['grid_z'][:, :])

        time1 = ma.getdata(ds['time'][:])
        velx = ma.getdata(velx[:, 0:, 0:])
        vely = ma.getdata(vely[:, 0:, 0:])

        if i == 0:
            N_t, N_files = len(velx), len(all_files)
            Y_m = np.zeros(N_t*N_files)
            time = np.zeros(N_t*N_files)

        Y_m1 = Average_Y(velx, vely, grid_x, grid_y)
        transition_instants, state_transition_count = switch_counts(
            Y_m1, time1, y_cut_plus, y_cut_minus)

        if state_transition_count != 0:
            print("Switch for run ", i+1,
                  " at time instants : ", transition_instants)
            total_transitions = total_transitions + state_transition_count

            for transition_instant in transition_instants:
                all_transition_instants = all_transition_instants + \
                    [time[i*N_t - 1] + transition_instant]

        Y_m[i*N_t:(i+1)*N_t] = Y_m1
        time[i*N_t:(i+1)*N_t] = time[i*N_t - 1] + time1

        i = i + 1

    print("Total count of switches: ", total_transitions)

    return Y_m, time, total_transitions, all_transition_instants


def switch_counts(Y_m, time, y_cut_plus, y_cut_minus, series_count=1, sample_N=0):
    """
    If series_count is specified more than 1, it needs sample size of the series
    , sample_N
    """

    #Initial states and counters, + is left, - is right
    if Y_m[0] > y_cut_plus:
        curr_state = 'left'
    elif Y_m[0] < y_cut_minus:
        curr_state = 'right'
    else:
        curr_state = 'center'

    state_transition_count = 0
    waiting_times = []
    counter = -1
    run = 0

    for t, y in zip(time, Y_m):

        counter = counter + 1
        if y > y_cut_plus:
            curr_state = 'left'
        elif y < y_cut_minus:
            curr_state = 'right'
        else:
            curr_state = 'center'

        if counter == 0:
            prev_state = curr_state
            prev_prev_state = curr_state
            transition_instants = []

        if curr_state != prev_state:
            switch = 1
        else:
            switch = 0

        #Check whether state goes from left to right or vice-versa
        if switch == 1:
            if curr_state != 'center' and curr_state != prev_prev_state and prev_prev_state != 'center':
                state_transition_count = state_transition_count + 1
                transition_instants.append(t)

            prev_prev_state = prev_state
        #Compute waiting time at the end of one continous series
        if counter == sample_N - 1:
            waiting_times = waiting_times + [transition_instants[0] - time[2000]*run] + (np.diff(transition_instants).tolist())
            #print(waiting_times)
            run = run + 1
                #[transition_instants[0]] + \
            counter = -1

        prev_state = curr_state

    return transition_instants, state_transition_count, waiting_times


def Plot_Ym(datadir, time, Y_m, all_transition_instants, y_cut_plus, y_cut_minus):

    ##Time series of combined runs
    fg_Y = plt.figure(figsize=(14, 6))
    font = {'family': 'normal', 'weight': 'normal', 'size': 32}
    plt.rc('font', **font)

    ax1 = fg_Y.add_subplot(111)
    ax1.plot(time, Y_m)
    ax1.set_xlabel('time (s)')
    ax1.set_ylabel('$Y_m$')
    plt.axhline(y=y_cut_plus, ls='--', c='g', label='state cut-off')
    plt.axhline(y=y_cut_minus, ls='--', c='g')

    for transition_instant in all_transition_instants[:-1]:
        plt.axvline(x=transition_instant, ls='--', c='r')

    #Use last instant for legend
    plt.axvline(x=all_transition_instants[-1],
                ls='--', c='r', label='switch instant')

    plt.legend(loc='lower right', bbox_to_anchor=(1, 1), ncol=3)
    plt.tight_layout()
    plt.savefig(datadir + "Ym_series.jpg", dpi=150)
    #plt.show()

#Plot drift, diffusion and pdf distribution


def Plot_DDP(datadir, N_bins, Y_m, stride, f_sample):

    ### Compute K-M coeffs + Edges ###
    dt = 1/f_sample

    bins = np.array([N_bins])

    ## The size is only valid for 1-D case
    KMc_exp = np.zeros((2, N_bins))

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
    ax1.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))

    ax2.plot(X_values, a_KM, 'ko')
    ax2.set_xlabel('Y')
    ax2.set_ylabel('a(Y)')
    ax2.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))

    ax3.plot(X_values, p_hist, 'ko')
    ax3.set_xlabel('Y')
    ax3.set_ylabel('p(Y)')

    plt.tight_layout()
    fg_dd_fit.subplots_adjust(wspace=0.2, hspace=0.4, top=0.85)

    fg_dd_fit.suptitle("Drift, diffusion and PDF plots")

    plt.savefig(datadir + "DDP.jpg", dpi=150)

    #plt.show()

def Lang_sim(f_a, g_a, X, *args, dt = 0.001, T=500, seed = 123, init = 'random'):

    #White noise parameters for noise term
    mu = 0.0
    std = 1.0
    seed = seed

    np.random.seed(seed)
    if init == 'random':
        y0 = np.random.choice(X[:,0])
    else :
        y0 = init

    #Other possible initalizations
    #y0 = np.mean(Ym_tot) + 1.5*np.std(Ym_tot)
    #y0 = Ym_tot[0,0]

    L = int(T/dt)

    y = np.empty(L, dtype=type(y0))
    t = np.empty(L, dtype=type(dt))
    y[0] = y0
    t[0] = 0

    #f_a = KMc_reg.get_exp_value(0)
    #g_a = KMc_reg.get_exp_value(1)


    #f_a = KM_temp.get_drift_fun(Xi=Xi[:,6])[0]
    #g_a = KM_temp.get_diff_fun(Xi=Xi[:,6])[0]

    count_up = 0
    count_down = 0
    for i in range(L-1):
        if y[i] > X[-1,0] :
            #print("Y is out of range:", y[i])
            count_up = count_up + 1
            y_index = -1
            f = f_a[y_index]
            g = g_a[y_index]

        elif y[i] < X[0,0] :
            #print("Y is out of range:", y[i])
            count_down = count_down + 1
            y_index = 0
            f = f_a[y_index]
            g = g_a[y_index]

        else:
            y_index = np.where(y[i] > X[:,0])[0][-1]
            f = ((f_a[y_index]) + (f_a[y_index+1]))/2
            g = (g_a[y_index] + g_a[y_index+1])/2

        dw = np.sqrt(dt)*np.random.normal(mu, std, 1)
        y[i+1] = y[i] + f*dt + g*dw
        t[i+1] = t[i] + dt

    print("Counts out of range up and down:", count_up, count_down)

    #Save simulation results at provided directory
    if len(args) != 0 :
        np.save(args[0] + "langevin_sim_spars_" + str(args[1]) + "_dt_" + str(dt) + "_T_" + str(T) + "_seed_" + str(seed), np.vstack((t, y)).T)

    return y, t

#Code for Markov test
def markov_test(Ym_tot, bins_hist):
    Ym_lagged = np.empty((Ym_tot.shape[0] - 2,3), dtype = Ym_tot.dtype)

    Ym_lagged[:,0] = Ym_tot[2:,0]
    Ym_lagged[:,1] = Ym_tot[1:-1,0]
    Ym_lagged[:,2] = Ym_tot[:-2,0]

    p1_joint_t3_t2_t1, bins_hist_3D = np.histogramdd(Ym_lagged, [bins_hist[0],bins_hist[0], bins_hist[0]], density=True)

    p_joint_t2_t1, bins_hist_2D = np.histogramdd(Ym_lagged[:,1:], [bins_hist[0],bins_hist[0]], density=True)

    p_joint_t3_t2, bins_hist_2D = np.histogramdd(Ym_lagged[:,:2], [bins_hist[0],bins_hist[0]], density=True)
    p_condn_t3_t2 = np.empty_like(p_joint_t3_t2)

    for i in range(len(p_joint_t3_t2)):
        p_condn_t3_t2[i,:] = p_joint_t3_t2[i,:]/p_hist

    p2_joint_t3_t2_t1 = np.einsum('ij,jk->ijk',p_condn_t3_t2 , p_joint_t2_t1)

    dx=bins_hist[0][1] - bins_hist[0][0]
    return utils.kl_divergence(p1_joint_t3_t2_t1, p2_joint_t3_t2_t1, dx=[dx, dx, dx], tol=1e-6)

if __name__ == "__main__":
    main()
