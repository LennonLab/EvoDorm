from __future__ import division
import FlowCytometryTools
from FlowCytometryTools import test_data_dir, test_data_file, FCMeasurement, FCPlate, ThresholdGate, PolyGate
import os, re, datetime, signal, collections
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.grid_search import GridSearchCV
from sklearn.neighbors import KernelDensity
from operator import itemgetter
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.animation as animation


data_path = os.path.expanduser('~/Box Sync/Dormancy_Visualization/data/all/')
git_path = os.path.realpath(__file__).rsplit('/', 2)[0]

samples = ['012', '015', '018', '021', '024', '036', '048', '060', '072', \
            '084', '096', '108', '120', '144', '168', '192', '216', '240', \
            '243', '246', '261', '264', '276', '288']

#samples = ['012']

class TimeoutException(Exception):   # Custom exception class
    pass

def timeout_handler(signum, frame):   # Custom signal handler
    raise TimeoutException


def makeFigFolders(analysis = 'RSG_DAPI'):
    if analysis == 'RSG':
        path = git_path + '/figures/flow_figs/RSG_DAPI/'
        if not os.path.exists(path):
            os.makedirs(path)
        for sample in samples:
            path_sample = path + 'Sample_' + sample + '/'
            if not os.path.exists(path_sample):
                os.makedirs(path_sample)
    #else:
    #    for sample in samples:

    #        rep_path = git_path + 'figures/flow_figs/' + analysis + '/' + strain + '/' + transfer + '/' + rep + '/'
    #        if not os.path.exists(rep_path):
    #            os.makedirs(rep_path)

def getDAPIgate(plate):
    As = ['A3', 'A4']
    cutoffs = []
    for A in As:
        DAPI = plate[A].data[['Pacific Blue-A']].values
        DAPI_gate = np.mean(DAPI) + (2*np.std(DAPI))
        cutoffs.append(DAPI_gate)
    cutoff = np.mean(cutoffs)
    return cutoff

def get_DAPI_RGF_Figs():
    As = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8']
    for sample in samples:
        print sample
        path = data_path + 'Sample_' + sample + '/'
        plate = FCPlate.from_dir(ID='Demo Plate', path = path, parser='name')
        plate = plate.dropna()
        plate = plate.transform('hlog', channels=['FSC-A', 'SSC-A', \
            'FSC PMT-A','PI (YG)-A', 'FITC-A', 'Pacific Blue-A', 'APC-A'])
        threshold = getDAPIgate(plate)
        for A in As:
            fig = plt.figure()
            DAPI = plate[A].data[['Pacific Blue-A']].values
            RSG = plate[A].data[['FITC-A']].values
            plt.scatter(DAPI, RSG)
            plt.axvline(x = threshold)
            flow_fig_path = git_path + '/figures/flow_figs/RSG_DAPI/Sample_' + sample + '/'
            if not os.path.exists(flow_fig_path):
                os.makedirs(flow_fig_path)
            name =  'RSG_DAPI_Sample_' + sample + '_' + A
            plt.savefig(flow_fig_path + name + '.png')

def get_dist_data():
    makeFigFolders(analysis = 'RSG')
    OUT1 = open(git_path +'/data/FlowSummary.txt', 'w+')
    print>> OUT1, 'Sample', 'Mean', 'Std', 'Median', 'Skew'
    for sample in samples:
        sample_name = 'Sample_' + sample
        path = data_path + 'Sample_' + sample + '/'
        files = os.listdir(path)
        remove = ['misc', '.DS_Store']
        files1 = [i for i in files if i not in remove]
        files_split = [re.split(r'[._]',x)[3] for x in files1]

        path = data_path + 'Sample_' + sample + '/'
        plate = FCPlate.from_dir(ID='Demo Plate', path = path, parser='name')
        plate = plate.dropna()
        plate = plate.transform('hlog', channels=['FSC-A', 'SSC-A', \
            'FSC PMT-A','PI (YG)-A', 'FITC-A', 'Pacific Blue-A', 'APC-A'])
        DAPI_threshold = getDAPIgate(plate)
        DAPI_gate = ThresholdGate(DAPI_threshold, 'Pacific Blue-A', region='above')
        intersection_gate = DAPI_gate
        DAPI_gated_A9 = plate['A8'].gate(intersection_gate)
        # RSG
        RSG = DAPI_gated_A9.data[['FITC-A']].values

        print>> OUT1, sample_name, str(np.mean(RSG)), str(np.std(RSG)), \
            str(np.median(RSG)), str(stats.skew(RSG)[0])

    OUT1.close()


def testCounts():
    tubes = ['A8']
    for sample in samples:
        print sample
        path = data_path + 'Sample_' + sample + '/'
        files = os.listdir(path)
        remove = ['misc', '.DS_Store']
        files1 = [i for i in files if i not in remove]
        files_split = [re.split(r'[._]',x)[3] for x in files1]
        if 'H' in files_split:
            for tube in tubes:
                file_folder_path = git_path + '/figures/flow_figs/' + tube
                if not os.path.exists(file_folder_path):
                    os.makedirs(file_folder_path)

                file_folder_sample_path = git_path + '/figures/flow_figs/' + \
                    tube +'/Sample_' + sample + '/'
                if not os.path.exists(file_folder_sample_path):
                    os.makedirs(file_folder_sample_path)

                plate = FCPlate.from_dir(ID='Demo Plate', path = path, parser='name')
                plate = plate.dropna()
                plate = plate.transform('hlog', channels=['FSC-A', 'SSC-A', \
                    'FSC PMT-A','PI (YG)-A', 'FITC-A', 'Pacific Blue-A', 'APC-A'])
                threshold = getDAPIgate(plate)
                gate = ThresholdGate(threshold, 'Pacific Blue-A', region='above')
                gated_sample = plate[tube].gate(gate)
                FSC_A = gated_sample.data[['FSC-A']].values
                SSC_A = gated_sample.data[['SSC-A']].values
                PI = gated_sample.data[['PI (YG)-A']].values
                FITC = gated_sample.data[['FITC-A']].values
                FSC_H = gated_sample.data[['FSC-H']].values
                SSC_H = gated_sample.data[['SSC-H']].values
                FSC_PMT_A = gated_sample.data[['FSC PMT-A']].values
                FSC_PMT_H = gated_sample.data[['FSC PMT-H']].values
                Pacific_blue = gated_sample.data[['Pacific Blue-A']].values
                APC = gated_sample.data[['APC-A']].values

                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                ax.scatter(FITC, Pacific_blue, APC , zdir='z', alpha = 0.2)
                ax.set_xlabel('FITC')
                ax.set_ylabel('Pacific blue')
                ax.set_zlabel('APC')

                name = tube + '_FITC_APC_PacificBlue_' + sample
                plt.savefig(file_folder_sample_path + name + '.png')
                plt.close()

                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.scatter(APC, FITC, alpha = 0.2)
                ax.set_xlabel('APC')
                ax.set_ylabel('FITC')

                name = tube +'_APC_FITC_' + sample
                plt.savefig(file_folder_sample_path + name + '.png')
                plt.close()


                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.scatter(Pacific_blue, FITC, alpha = 0.2)
                ax.set_xlabel('Pacific blue')
                ax.set_ylabel('FITC')

                name = tube + '_PacificBlue_FITC_' + sample
                plt.savefig(file_folder_sample_path + name + '.png')
                plt.close()

                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.scatter(FSC_A, SSC_A, alpha = 0.2)
                ax.set_xlabel('FSC-A')
                ax.set_ylabel('SSC-A')

                name = tube + '_FSC_A_SSC_A_' + sample
                plt.savefig(file_folder_sample_path + name + '.png')
                plt.close()

                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.scatter(FSC_H, SSC_H, alpha = 0.2)
                ax.set_xlabel('FSC-H')
                ax.set_ylabel('SSC-H')
                #plt.axhline(y = 250000, linewidth=2, color='darkgrey', ls='--')

                name = tube + '_FSC_H_SSC_H_' + sample
                plt.savefig(file_folder_sample_path + name + '.png')
                plt.close()

                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.scatter(FSC_PMT_A, FSC_PMT_H, alpha = 0.2)
                ax.set_xlabel('FSC PMT-A')
                ax.set_ylabel('FSC PMT-H')

                name = tube + '_FSC_PMT_A_FSC_PMT_H_' + sample
                plt.savefig(file_folder_sample_path + name + '.png')
                plt.close()


def get_RSG_hist():
    folder_path = git_path + '/figures/flow_figs/RSG'
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    fig1 = plt.figure()
    for sample in samples:
        print sample
        path = data_path + 'Sample_' + sample + '/'
        plate = FCPlate.from_dir(ID='Demo Plate', path = path, parser='name')
        plate = plate.dropna()
        plate = plate.transform('hlog', channels=['FSC-A', 'SSC-A', \
            'FSC PMT-A','PI (YG)-A', 'FITC-A', 'Pacific Blue-A', 'APC-A'])
        threshold = getDAPIgate(plate)
        gate = ThresholdGate(threshold, 'Pacific Blue-A', region='above')
        gated_sample = plate['A8'].gate(gate)
        RSG = gated_sample.data[['FITC-A']].values

        fig, ax = plt.subplots()
        ax.hist(RSG, 30, fc='gray', histtype='stepfilled', alpha=0.3, normed=True)
        axes.set_xlim([0,10000])
        plt.savefig(git_path + '/figures/flow_figs/RSG/RSG_Sample_' + sample + '.png')
        plt.close()



def get_RSG_movie():
    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    number_of_frames = len(samples)
    for sample in samples:
        print sample
        path = data_path + 'Sample_' + sample + '/'
        plate = FCPlate.from_dir(ID='Demo Plate', path = path, parser='name')
        plate = plate.dropna()
        plate = plate.transform('hlog', channels=['FSC-A', 'SSC-A', \
            'FSC PMT-A','PI (YG)-A', 'FITC-A', 'Pacific Blue-A', 'APC-A'])
        threshold = getDAPIgate(plate)
        gate = ThresholdGate(threshold, 'Pacific Blue-A', region='above')
        gated_sample = plate['A8'].gate(gate)
        RSG = gated_sample.data[['FITC-A']].values


        l, = plt.plot([], [], 'r-')
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.xlabel('x')
        plt.title('test')
        line_ani = animation.FuncAnimation(fig1, update_line, 25, fargs=(data, l),
                                           interval=50, blit=True)
        line_ani.save('lines.mp4', writer=writer)


        def update_hist(num, data):
            plt.cla()
            plt.hist(data[num])

        fig = plt.figure()
        hist = plt.hist(RSG, 30, fc='gray', histtype='stepfilled', alpha=0.3, normed=True)

        #animation1 = animation.FuncAnimation(fig, update_hist, number_of_frames, fargs=(data, ) )
        #plt.savefig(git_path + '/figures/flow_figs/RSG/RSG_Sample_' + sample + '.png')
        #plt.close()

def movie_test():
    #print animation.writers
    #Writer = animation.writers['ffmpeg']
    #writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    n = 100
    number_of_frames = 100
    data = np.random.rand(n, number_of_frames)

    def update_hist(num, data):
        plt.cla()
        plt.hist(data[num])

    fig = plt.figure()
    hist = plt.hist(data[0])
    print data[0]

    animation1 = animation.FuncAnimation(fig, update_hist, frames= number_of_frames, fargs=(data, ) )
    movie_name = git_path + '/figures/lines.mp4'
    animation1.save(movie_name, fps=4)

#get_dist_data()
#testCounts()
#makeFigFolders()
#get_DAPI_RGF_Figs()
#get_KDE_RSG()

movie_test()
