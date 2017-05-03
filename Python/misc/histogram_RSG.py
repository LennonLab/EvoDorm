"""
==================
Animated histogram
==================

This example shows how to use a path patch to draw a bunch of
rectangles for an animated histogram.

"""
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path
import matplotlib.animation as animation
from FlowCytometryTools import test_data_dir, test_data_file, FCMeasurement, FCPlate, ThresholdGate, PolyGate


samples = ['003', '006', '009', '012', '015', '018', '021', '024', '036', '048', '060', '072', \
            '084', '096', '108', '120', '144', '168', '192', '216', '240', \
            '243', '246', '261', '264', '276', '288', '300', '312', '324', '336', \
            '360', '384', '408', '432', '456']

git_path = os.path.realpath(__file__).rsplit('/', 2)[0]
data_path = os.path.expanduser('~/Box Sync/Dormancy_Visualization/data/all/')

def getDAPIgate(plate):
    As = ['A3', 'A4']
    cutoffs = []
    for A in As:
        DAPI = plate[A].data[['Pacific Blue-A']].values
        DAPI_gate = np.mean(DAPI) + (2*np.std(DAPI))
        cutoffs.append(DAPI_gate)
    cutoff = np.mean(cutoffs)
    return cutoff


fig, ax = plt.subplots()

# histogram our data with numpy

path1 = data_path + 'Sample_' + samples[0] + '/'
plate1 = FCPlate.from_dir(ID='Demo Plate', path = path1, parser='name')
plate1 = plate1.dropna()
plate1 = plate1.transform('hlog', channels=['FSC-A', 'SSC-A', \
    'FSC PMT-A','PI (YG)-A', 'FITC-A', 'Pacific Blue-A', 'APC-A'])
threshold = getDAPIgate(plate1)
gate1 = ThresholdGate(threshold, 'Pacific Blue-A', region='above')
gated_sample1 = plate1['A8'].gate(gate1)
RSG1 = gated_sample1.data[['FITC-A']].values
n, bins = np.histogram(RSG1, 40, normed=True)

# get the corners of the rectangles for the histogram
left = np.array(bins[:-1])
right = np.array(bins[1:])
bottom = np.zeros(len(left))
top = bottom + n
nrects = len(left)

# here comes the tricky part -- we have to set up the vertex and path
# codes arrays using moveto, lineto and closepoly

# for each rect: 1 for the MOVETO, 3 for the LINETO, 1 for the
# CLOSEPOLY; the vert for the closepoly is ignored but we still need
# it to keep the codes aligned with the vertices
nverts = nrects*(1 + 3 + 1)
verts = np.zeros((nverts, 2))
codes = np.ones(nverts, int) * path.Path.LINETO
codes[0::5] = path.Path.MOVETO
codes[4::5] = path.Path.CLOSEPOLY
verts[0::5, 0] = left
verts[0::5, 1] = bottom
verts[1::5, 0] = left
verts[1::5, 1] = top
verts[2::5, 0] = right
verts[2::5, 1] = top
verts[3::5, 0] = right
verts[3::5, 1] = bottom

barpath = path.Path(verts, codes)
patch = patches.PathPatch(
    barpath, facecolor='green', edgecolor='yellow', alpha=0.5)
ax.add_patch(patch)

ax.set_xlim(left[0], right[-1])
ax.set_ylim(0, 0.001)
time_text = ax.text(0.85, 0.9, '', transform=ax.transAxes)
ax.set_title('The distribution of reductase \n acivity in a microbial population')
ax.set_xlabel('Reductase-associated fluorescence')


def animate(i):
    sample = samples[i]
    path = data_path + 'Sample_' + sample + '/'
    plate = FCPlate.from_dir(ID='Demo Plate', path = path, parser='name')
    plate = plate.dropna()
    plate = plate.transform('hlog', channels=['FSC-A', 'SSC-A', \
        'FSC PMT-A','PI (YG)-A', 'FITC-A', 'Pacific Blue-A', 'APC-A'])
    threshold = getDAPIgate(plate)
    gate = ThresholdGate(threshold, 'Pacific Blue-A', region='above')
    gated_sample = plate['A8'].gate(gate)
    RSG = gated_sample.data[['FITC-A']].values

    # simulate new data coming in
    n, bins = np.histogram(RSG, 40, normed=True)
    top = bottom + n
    verts[1::5, 1] = top
    verts[2::5, 1] = top
    #t = ax.annotate(sample + ' hours',(7900,1200)) # add text
    time_text.set_text(sample + ' hours' )
    return [patch, time_text]

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ani = animation.FuncAnimation(fig, animate, len(samples), repeat=False, blit=True)
movie_name = git_path + '/figures/hist.gif'
ani.save(movie_name, fps=2, writer='imagemagick')
