from __future__ import division
import pandas as pd
import numpy as np
import  matplotlib.pyplot as plt
import os, itertools, math
import numpy.polynomial.polynomial as poly
import matplotlib.lines as mlines
from scipy.stats import gaussian_kde
from scipy import stats


mydir = os.path.expanduser("~/GitHub/EvoDormReview/")


def dist_core_genome_clean():
    OUT = open(mydir + 'data/G50_U0.001_N200.txt', 'w+')
    print>> OUT, 'r', 'pan_genome_mean', 'core_genome_mean', 'pan_genome_std', \
        'core_genome_std'
    for r_file in os.listdir(mydir + 'data/r_range'):
        if r_file.endswith(".txt"):
            IN = pd.read_csv(mydir + 'data/r_range/' + r_file, sep = ' ')
            pan = IN['pan_genome'].values
            core = IN['core_genome'].values
            total = pan + core
            f_pan = pan / total
            f_core = core / total
            f_pan_mean = np.mean(f_pan)
            f_core_mean = np.mean(f_core)
            f_pan_std = np.std(f_pan)
            f_core_std = np.std(f_core)
            r = float(r_file.split('_')[2][1:])
            print>> OUT, r, f_pan_mean, f_core_mean, f_pan_std, f_core_std
    OUT.close()

def dist_core_genome_fig():
    IN = pd.read_csv(mydir + 'data/G50_U0.001_N200.txt', sep = ' ')
    sort = IN.sort('r', ascending=False)
    r = np.log10(sort['r'].values)
    fig, ax = plt.subplots()
    colors = ['#87CEEB',  '#FF6347']
    plt.grid()
    ax.plot(r, sort['pan_genome_mean'],  lw = 2, color = colors[0], label = 'Distributed genome')
    ax.fill_between(r, sort['pan_genome_mean']+sort['pan_genome_std'], \
        sort['pan_genome_mean']-sort['pan_genome_std'], facecolor=colors[0], \
        alpha=0.5)
    ax.plot(r, sort['core_genome_mean'],  lw = 2, color = colors[1], lettabel = 'Core genome')
    ax.fill_between(r, sort['core_genome_mean']+sort['core_genome_std'], \
        sort['core_genome_mean']-sort['core_genome_std'], facecolor=colors[1], \
        alpha=0.5)
    ax.set_xlabel('Recombination rate, ' + r'$log_{10}$', fontsize=20)
    ax.set_ylabel( 'Frequency of genes', fontsize=20)
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, \
           ncol=2, mode="expand", borderaxespad=0.)
    fig.savefig(mydir + 'figures/test.png', \
    bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.
    """
    if not exponent:
        exponent = int(np.floor(np.log10(abs(num))))
    coeff = round(num / float(10**exponent), decimal_digits)
    if not precision:
        precision = decimal_digits

    return r"{0:.{2}f}\cdot10^{{{1:d}}}".format(coeff, exponent, precision)


def Fig2(theta2 = 0, seq_length = 100):
    fig = plt.figure()
    for i in range(2):
        ax = fig.add_subplot(2, 1, i+1)
        if i == 0:
            N = 1000
            theta1s = [2* N* 0.0000005, 2* N*0.00000005, 2* N*0.000000005]
            colors = ['#FF6347', '#FFA500', '#87CEEB']
            K = np.logspace(-3, 3, num = 1000, base=10.0)
            M = [(N/k) for k in K]
            for i, theta1 in enumerate(theta1s):
                y = []
                for K_i in K:
                    term1 = theta1 + (theta1 / K_i)
                    term2 = (1 + (1 / K_i))  * (theta2 / K_i)
                    pi = term1 + term2
                    y.append(pi)
                theta1_SN = sci_notation(theta1)
                ax.plot(M, np.asarray(y), 'k-', color = colors[i], label= r'$\theta = {{{}}}$'.format(theta1_SN), linewidth=3)
            ax.set_ylim([0.000001, 1])
            ax.set_xscale('log', basex=10)
            ax.set_yscale('log', basey=10)

            plt.axvline(x = 1000, c = 'grey', linestyle = '--', lw = 3)
            ax.legend(loc='upper left', prop={'size':12})
            ax.set_xlabel('Avg. time in seed bank (generations), '+ r'$log_{10}$' , fontsize=18)
            ax.set_ylabel('Nucleotide diversity (' + r'$\pi$' + '), ' + r'$log_{10}$', \
                fontsize=14, labelpad= 16)

        elif i == 1:
            Ms = [10, 100, 1000]
            Ex = [2.02, 2.2, 4]
            colors = ['#87CEEB', '#FFA500', '#FF6347']
            x = np.arange(1,10001)
            for i, M in enumerate(Ms):
                IN = pd.read_csv(mydir + 'data/Fig2b/merged/G10000_S100_N1000_M' + str(M) + '_c10_Pi.txt', \
                    sep = ' ')
                IN_mean = IN.mean(axis = 1).values / seq_length
                IN_std = IN.std(axis = 1)
                ax.plot(x, IN_mean, lw=2, color=colors[i], \
                    label = 'Time in seed bank = ' + r'${{{}}}$'.format(M) ,alpha = 0.9)
                ax.axhline(y = Ex[i] / seq_length, color=colors[i], ls = '--', lw = 3)
            ax.legend(loc='upper left', prop={'size':12})
            ax.set_xlabel('Time (generations)', fontsize=18)
            ax.set_ylabel('Nucleotide diversity (' + r'$\pi$' + ')', fontsize=14, \
                labelpad= 16)
            ax.set_ylim([0,0.05])
    #fig.text(.05, 0.5, 'Expected nucleotide diversity (' + r'$\pi$' + ')' + r'$log_{10}$', \
    #    fontsize=18, rotation='vertical', \
    #horizontalalignment='center', verticalalignment='center')
    fig.text(0.15, 0.95, 'a)',  fontsize=14,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.15, 0.475, 'b)',  fontsize=14,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.tight_layout()
    fig.savefig(mydir + 'figures/Fig2.png', \
            bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def Fig3(N = 1000, M = 10000, s = 0.1):
    fig = plt.figure()
    for i in range(2):
        if i == 0:
            ax = fig.add_subplot(2, 1, i + 1)
            Cs = [1, 10, 100]
            maxg = 0
            colors = ['#FF6347', '#FFA500', '#87CEEB']
            for j, c in enumerate(Cs):
                g = []
                p = []
                df = pd.read_csv(mydir + 'data/Fig3a/sweep_N_' + str(N) + '_M_1000_c_' + \
                    str(c) + '_s_' + str(s) + '.txt', header = None)
                df = df.fillna(1.0)
                for index, row in df.iterrows():
                    g.extend(row[~row.isnull()].index.values)
                    p.append(row[~row.isnull()].values)
                #g_plot = [np.mean(x) for x in zip(*g)]
                #print itertools.izip_longest(*p)
                p_plot_mean = []
                p_plot_std = []
                for i in itertools.izip_longest(*p):
                    i = np.asarray(i)
                    i = i[i != np.array(None)]
                    p_plot_mean.append( np.mean(i))
                    p_plot_std.append( np.std(i))
                p_plot_mean = np.asarray(p_plot_mean)
                p_plot_std = np.asarray(p_plot_std)
                ax.plot(range(1, len(p_plot_mean )+1),p_plot_mean, label='Time in seed bank = ' + str(int(1 / (c/M))), lw = 2, color = colors[j])
                ax.fill_between(range(1, len(p_plot_mean )+1), p_plot_mean+p_plot_std, p_plot_mean-p_plot_std, facecolor=colors[j], alpha=0.5)

                if max(g) > maxg:
                    maxg = max(g)
            #ax.set_xlim([0, maxg])
            ax.set_ylim([0, 1])
            ax.legend(loc='upper left', fontsize = 12)
            ax.set_xscale('log', basex=10)
            plt.axhline(y = 0.5, c = 'grey', linestyle = '--', lw = 3)
            ax.set_xlabel('Time (generations), ' + r'$log_{10}$', fontsize=20)
            ax.set_ylabel('Frequency of favored allele', fontsize=14)

        elif i == 1:
            ax = fig.add_subplot(2, 1, i+1)
            colors = ['#FF6347', '#FFA500', '#87CEEB']
            pop_type = {'N': 'Active', 'M': 'Dormant'}
            IN_N = pd.read_csv(mydir + 'data/Fig3b/T_fix_N_sweep_N_1000_M_10000_s_0.1_r_100.txt', sep = ' ')
            IN_M = pd.read_csv(mydir + 'data/Fig3b/T_fix_M_sweep_N_1000_M_10000_s_0.1_r_100.txt', sep = ' ')
            df_add = IN_N.add(IN_M, fill_value=0)
            df = df_add.divide(2, axis=0)
            Cs = df.columns.values.astype(float)
            timeInSb = 1 / (Cs / M)
            timeInSb_N = 1 / (IN_N.columns.values.astype(float) / M)
            timeInSb_M = 1 / (IN_M.columns.values.astype(float) / M)
            means = []
            std = []
            means_N = []
            std_N = []
            means_M = []
            std_M = []

            for column in df:
                data = df[column].values.astype(float)
                means.append(np.mean(data))
                std.append(np.std(data))
            for column in IN_N:
                data = IN_N[column].values.astype(float)
                means_N.append(np.mean(data))
                std_N.append(np.std(data))
            for column in IN_M:
                data = IN_M[column].values.astype(float)
                means_M.append(np.mean(data))
                std_M.append(np.std(data))

            means = np.asarray(means)
            std = np.asarray(std)
            means_N = np.asarray(means_N)
            std_N = np.asarray(std_N)
            means_M = np.asarray(means_M)
            std_M = np.asarray(std_M)

            ax.plot(timeInSb, means,  lw = 2, label='Active and dormant', \
                color = colors[0])
            ax.fill_between(timeInSb, means+std, means-std, facecolor=colors[0], alpha=0.5)

            ax.plot(timeInSb_N, means_N,  lw = 2, label='Active', \
                color = colors[1])
            ax.fill_between(timeInSb_N, means_N+std_N, means_N-std_N, facecolor=colors[1], alpha=0.5)

            ax.plot(timeInSb_M, means_M,  lw = 2, label='Dormant', \
                color = colors[2])
            ax.fill_between(timeInSb_M, means_M+std_M, means_M-std_M, facecolor=colors[2], alpha=0.5)

            ax.set_xlabel('Average time in seed bank, ' + r'$log_{10}$', fontsize=20)
            ax.legend(loc='upper left', fontsize = 12)
            ax.set_ylabel( r'$T_{fix}$', fontsize=20)
            plt.axvline(x = N, c = 'grey', linestyle = '--', lw = 3)
            plt.axhline(y = 117.5107, c = 'blue', linestyle = '--', lw = 3)
            ax.set_xscale('log', basex=10)
            ax.set_yscale('log', basey=10)
            ax.set_ylim([10, 100000])

    fig.text(0.14, 0.955, 'a)',  fontsize=14,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.14, 0.48, 'b)',  fontsize=14,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.tight_layout()
    fig.savefig(mydir + 'figures/Fig3.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def Fig4(subpop = 'all'):
    fig = plt.figure()
    for i in range(2):
        if i == 0:
            ax = fig.add_subplot(2, 1, i + 1)
            IN = pd.read_excel(mydir + 'data/Fig4/evo12597-sup-0004-Table-S2.xlsx')
            IN.columns = ['Taxon', 'NCBI', 'SporeGenes', 'dS', 'BranchLength', 'CodonBias', 'SporeForming']
            x = IN['SporeGenes'].values
            x = np.log10(x)
            y = IN['BranchLength'].values
            SporeForming = IN['SporeForming'].values
            colors = [ '#FF8C00' if i == 'N' else '#4169E1' for i in SporeForming ]
            NS_N = zip(colors, x, y)
            NS = [z[1:] for z in NS_N if z[0] == '#FF8C00']
            NS_x = [z[0] for z in NS]
            NS_y = [z[1] for z in NS]
            S  = [z[1:] for z in NS_N if z[0] == '#4169E1' ]
            S_x = [z[0] for z in S]
            S_y = [z[1] for z in S]
            ax.scatter(NS_x, NS_y, c='#FF8C00', marker='o', label='Non-spore forming')
            ax.scatter(S_x, S_y, c='#4169E1', marker='o', label='Spore forming')
            ax.legend(loc='upper right')
            slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
            predict_y = intercept + slope * x
            pred_error = y - predict_y
            degrees_of_freedom = len(x) - 2
            residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
            ax.plot(x, predict_y, 'k-', lw = 2)

            ax.set_xlabel('Number of sporulation genes, ' + r'$log_{10}$', fontsize=18)
            ax.set_ylabel('Evolutionary distance \n (root-to-tip distance)', fontsize=16)

        elif i == 1:
            ax = fig.add_subplot(2, 1, i + 1)
            cs = np.linspace(1, 100, num = 50, endpoint=True)
            cs = np.rint(cs)
            distances = []
            cs_x = []
            distances_c_x = []
            for c in cs:
                c = int(c)
                IN_path = mydir + 'data/Fig4/Fig4_sim/Fig4_sim_c_' + str(c) +'.txt'
                IN = pd.read_csv(IN_path, sep = ' ', header = 'infer')
                #data_dict[avg_time] = IN.evol_distance.values
                if subpop == 'all':
                    distance = IN.evol_distance.values
                elif subpop == 'N':
                    distance = IN.evol_distance.values
                elif subpop == 'M':
                    distance = IN.evol_distance.values
                distances.append(distance)
                distances_c_x.extend(distance)
                cs_x.extend(list(itertools.repeat(c, len(distance))))
            average_times = [ 1 / (c/100) for c in cs]
            average_times_cs_x = [ 1 / (c/100) for c in cs_x]
            average_times = np.log10(average_times)
            distances_mean = [np.mean(x) for x in distances]
            distances_std = [np.std(x) for x in distances]
            x =  np.log10(np.asarray(average_times_cs_x))
            y = np.asarray(distances_c_x)
            xy = np.vstack([x, y])
            z = gaussian_kde(xy)(xy)
            ax.scatter(x, y, c=z, s=100, edgecolor='', \
                cmap = 'viridis', alpha = 0.8)

            slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
            predict_y = intercept + slope * x
            pred_error = y - predict_y
            degrees_of_freedom = len(x) - 2
            residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
            ax.plot(x, predict_y, 'k-', lw = 2)
            ax.set_xlim([0, 2.05])
            ax.set_ylim([0, 0.5])
            ax.set_xlabel('Average time in seed bank (generations), ' \
                + r'$log_{10}$', fontsize=18)
            ax.set_ylabel('Evolutionary distance \n (JC69-corrected distance)', fontsize=16)

    if subpop == 'all':
        fig_name = mydir + 'figures/Fig4.png'
    elif subpop == 'N':
        fig_name = mydir + 'figures/Fig4_N.png'
    elif subpop == 'M':
        fig_name = mydir + 'figures/Fig4_M.png'
    fig.text(0.15, 0.97 , 'a)',  fontsize=14,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.15, 0.475, 'b)',  fontsize=14,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.tight_layout()
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()




#Fig4()
