from __future__ import division
import pandas as pd
import numpy as np
import  matplotlib.pyplot as plt
import os, itertools, math
import numpy.polynomial.polynomial as poly
import matplotlib.lines as mlines
from scipy.stats import gaussian_kde
from scipy import stats
from FlowCytometryTools import FCPlate, ThresholdGate


mydir = os.path.realpath(__file__).rsplit('/', 2)[0]


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
                IN = pd.read_csv(mydir + '/data/Fig2b/merged/G10000_S100_N1000_M' + str(M) + '_c10_Pi.txt', \
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
    fig.savefig(mydir + '/figures/Fig2.png', \
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
                df = pd.read_csv(mydir + '/data/Fig3a/sweep_N_' + str(N) + '_M_1000_c_' + \
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
            IN_N = pd.read_csv(mydir + '/data/Fig3b/T_fix_N_sweep_N_1000_M_10000_s_0.1_r_100.txt', sep = ' ')
            IN_M = pd.read_csv(mydir + '/data/Fig3b/T_fix_M_sweep_N_1000_M_10000_s_0.1_r_100.txt', sep = ' ')
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
                color = colors[1])
            ax.fill_between(timeInSb, means+std, means-std, facecolor=colors[1], alpha=0.5)

            ax.plot(timeInSb_N, means_N,  lw = 2, label='Active', \
                color = colors[2])
            ax.fill_between(timeInSb_N, means_N+std_N, means_N-std_N, facecolor=colors[2], alpha=0.5)

            ax.plot(timeInSb_M, means_M,  lw = 2, label='Dormant', \
                color = colors[0])
            ax.fill_between(timeInSb_M, means_M+std_M, means_M-std_M, facecolor=colors[0], alpha=0.5)

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
    fig.savefig(mydir + '/figures/Fig3.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


def Fig4(subpop = 'all'):
    fig = plt.figure()
    for i in range(2):
        if i == 0:
            ax = fig.add_subplot(2, 1, i + 1)
            IN = pd.read_excel(mydir + '/data/Fig4/evo12597-sup-0004-Table-S2.xlsx')
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
                IN_path = mydir + '/data/Fig4/Fig4_sim/Fig4_sim_c_' + str(c) +'.txt'
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
        fig_name = mydir + '/figures/Fig4.png'
    elif subpop == 'N':
        fig_name = mydir + '/figures/Fig4_N.png'
    elif subpop == 'M':
        fig_name = mydir + '/figures/Fig4_M.png'
    fig.text(0.15, 0.97 , 'a)',  fontsize=14,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.text(0.15, 0.475, 'b)',  fontsize=14,
        horizontalalignment='center',
        verticalalignment='center', fontweight='bold')
    fig.tight_layout()
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def getDAPIgate(plate):
    As = ['A3', 'A4']
    cutoffs = []
    for A in As:
        DAPI = plate[A].data[['Pacific Blue-A']].values
        DAPI_gate = np.mean(DAPI) + (2*np.std(DAPI))
        cutoffs.append(DAPI_gate)
    cutoff = np.mean(cutoffs)
    return cutoff

def Box1Fig():
    fig = plt.figure()

    path_006 = mydir + '/data/Box1Fig/Sample_006/'
    path_012 = mydir + '/data/Box1Fig/Sample_012/'
    path_264 = mydir + '/data/Box1Fig/Sample_264/'
    plate_006 = FCPlate.from_dir(ID='Demo Plate', path = path_006, parser='name')
    plate_012 = FCPlate.from_dir(ID='Demo Plate', path = path_012, parser='name')
    plate_264 = FCPlate.from_dir(ID='Demo Plate', path = path_264, parser='name')
    plate_006 = plate_006.dropna()
    plate_012 = plate_012.dropna()
    plate_264 = plate_264.dropna()
    plate_006 = plate_006.transform('hlog', channels=['FSC-A', 'SSC-A', \
        'FSC PMT-A','PI (YG)-A', 'FITC-A', 'Pacific Blue-A', 'APC-A'])
    plate_012 = plate_012.transform('hlog', channels=['FSC-A', 'SSC-A', \
        'FSC PMT-A','PI (YG)-A', 'FITC-A', 'Pacific Blue-A', 'APC-A'])
    plate_264 = plate_264.transform('hlog', channels=['FSC-A', 'SSC-A', \
        'FSC PMT-A','PI (YG)-A', 'FITC-A', 'Pacific Blue-A', 'APC-A'])
    threshold_006 = getDAPIgate(plate_006)
    threshold_012 = getDAPIgate(plate_012)
    threshold_264 = getDAPIgate(plate_264)
    gate_006 = ThresholdGate(threshold_006, 'Pacific Blue-A', region='above')
    gate_012 = ThresholdGate(threshold_012, 'Pacific Blue-A', region='above')
    gate_264 = ThresholdGate(threshold_264, 'Pacific Blue-A', region='above')
    gated_sample_006 = plate_006['A8'].gate(gate_006)
    gated_sample_012 = plate_012['A8'].gate(gate_012)
    gated_sample_264 = plate_264['A8'].gate(gate_264)
    RSG_006 = gated_sample_006.data[['FITC-A']].values
    RSG_012 = gated_sample_012.data[['FITC-A']].values
    RSG_264 = gated_sample_264.data[['FITC-A']].values

    #colors = ['#FF6347', '#FFA500', '#87CEEB']

    plt.hist(RSG_006, 40, fc='#87CEEB', histtype='bar', alpha=0.5, normed=True)
    plt.hist(RSG_012, 40, fc='#FFA500', histtype='bar', alpha=0.5, normed=True)
    plt.hist(RSG_264, 40, fc='#FF6347', histtype='bar', alpha=0.5, normed=True)

    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.title('The distribution of reductase \n acivity in a microbial population', \
        fontsize = 20, weight = 'heavy')
    plt.xlabel('Reductase-associated fluorescence', fontsize = 18)
    plt.ylabel('Probability', fontsize = 18)
    plt.arrow(2400, 0.00075, 3300, 0, width=0.00004, \
        head_width=0.00012, head_length=450,  length_includes_head=True, \
        shape='full', color = '#87CEEB')
    plt.text(2800, 0.0008, 'Resucitation', color = '#87CEEB', fontsize = 18, weight = 'heavy')
    plt.arrow(5250, 0.00064, -3300, 0, width=0.00004, \
        head_width=0.00012, head_length=450,  length_includes_head=True, \
        shape='full', color = '#FF6347')
    plt.text(2975, 0.00055 , 'Dormancy', color = '#FF6347', fontsize = 18, weight = 'heavy')

    plt.xlim(0, 8000)
    plt.ylim(0, 0.001)
    fig.tight_layout()
    fig_name = mydir + '/figures/Box1Fig.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


Box1Fig()
