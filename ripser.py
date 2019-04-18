'''
Created by Thiago de Melo
'''


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import subprocess
#import os

def circle_pts(N, R=1):
    theta_list = np.random.random_sample(N)
    pts = np.zeros((N,2))
    #print theta_list
    for i in range(len(theta_list)):
        pts[i,0] = R*np.cos(2*np.pi*theta_list[i])
        pts[i,1] = R*np.sin(2*np.pi*theta_list[i])
    return pts

def sphere_pts(N, R=1):
    angle_list = 2*np.pi * np.random.sample((N,2))
    #phi_list   = np.random.random_sample(N)
    pts = np.zeros((N,3))
    for i in range(len(angle_list)):
        pts[i,0] = R*np.cos(angle_list[i,0])*np.cos(angle_list[i,1])
        pts[i,1] = R*np.cos(angle_list[i,0])*np.sin(angle_list[i,1])
        pts[i,2] = R*np.sin(angle_list[i,0])
    return pts

def figure_eight_pts(N, a=1):
    theta_list = 2 * np.pi * np.random.sample(N)
    pts = np.zeros((N,2))
    print theta_list
    for i in range(len(theta_list)):
        pts[i,0] = a * np.cos(theta_list[i]) * np.sqrt(2*np.cos(2*theta_list[i]))
        pts[i,1] = a * np.sin(theta_list[i]) * np.sqrt(2*np.cos(2*theta_list[i]))
    return pts

def annulus_pts(N, R=2, r=1):
    theta_list = np.random.random_sample(N)
    radius_list = r + np.random.random_sample(N) * (R-r)
    pts = np.zeros((N,2))
    for i in range(len(theta_list)):
        pts[i,0] = radius_list[i] * np.cos(2*np.pi*theta_list[i])
        pts[i,1] = radius_list[i] * np.sin(2*np.pi*theta_list[i])
    return pts

def random_pts(N,d):
    pts = np.random.random_sample((N, d))
    return pts

def cube_pts(N):
    npts = N/6
    faces = {}
    for i in range(3):
        data0 = np.random.random((npts,3))
        data1 = np.random.random((npts,3))
        data0[:,i] = 0
        data1[:,i] = 1
        faces[i]   = data0
        faces[i+3] = data1
    cube = np.concatenate([faces[i] for i in range(6)])
    return cube

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 5))  # outward by 10 points
            spine.set_smart_bounds(True)
        else:
            spine.set_color('none')  # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])

def plot_bar(p,q,c='b',linestyle='-'):
    plt.plot([p[0],q[0]],[p[1],q[1]], c=c,linestyle=linestyle, linewidth=1)


''' DATA SET SECTION '''
''' cant use dash - in variable name (use _ instead)'''
VR = './examples/VR.point_cloud'
mds_weight = './examples/mds-weight.point_cloud'
mds_size = './examples/mds-size.point_cloud'
input_file = mds_weight
pts_file = './pts.point_cloud'

'''choose pts below'''
#pts = np.loadtxt(input_file, delimiter=',')
pts = circle_pts(50, 2)
#pts = annulus_pts(400, 4, 1)
#pts = random_pts(100,3)
#pts = cube_pts(30)
#pts = figure_eight_pts(400, 3)
#pts = sphere_pts(100)

np.savetxt(pts_file, pts, delimiter=' ', fmt='%.3f')
num_pts, dim_pts = pts.shape

#print pts

#plt.scatter(pts[:,0],pts[:,1])
#plt.show()
#quit()

''' RIPSER SECTION '''
ripser_bin = ['/home/thiago/Dropbox/Programacao/ripser/ripser']
ripser_input_file = pts_file
ripser_output_file = 'ripser.out'
ripser_format = 'point-cloud'
ripser_threshold = 1.4
ripser_dim = 2
ripser_opt = ['--format', ripser_format, '--threshold', str(ripser_threshold), '--dim', str(ripser_dim), '--output', ripser_output_file, ripser_input_file]
ripser_cmd = ripser_bin + ripser_opt

print "Executing Ripser..."
subprocess.call(ripser_cmd),
print "Done!"

#quit()

''' DIAGRAMS SECTION '''
cols_name = ['dim', 'birth', 'death']
df = pd.read_csv(ripser_output_file, delim_whitespace = True, header = None, names = cols_name)
dimensions = df.drop_duplicates('dim')['dim'].tolist()
birth_max = df['birth'].max()
death_max = df['death'].max()
infinity = ripser_threshold #1.05 * np.maximum(birth_max,death_max)
print "Creating persistent diagram for dimensions", dimensions

diagrams = {} 
diagrams_inf = {}
bar_num = {}
bar_inf_num = {}
for dim in dimensions:
    diagrams[dim] = df[ (df.dim == dim) & (df.death != -1) ]
    diagrams[dim] = diagrams[dim].sort_values(['birth','death'], ascending=[True,True]).reset_index()
    diagrams_inf[dim] = df[ (df.dim == dim) & (df.death == -1) ]
    diagrams_inf[dim] = diagrams_inf[dim].sort_values(['birth','death'], ascending=[True,True]).reset_index()
    bar_num[dim] = diagrams[dim].shape[0]
    bar_inf_num[dim] = diagrams_inf[dim].shape[0]
    #print dim, bar_num[dim], bar_inf_num[dim]
    #print diagrams[dim]


''' PLOTS SECTION '''

#'''
fig = plt.figure()

# dataset
if dim_pts >= 3:
    ax = fig.add_subplot(1,2,2, projection='3d')
    ax.plot3D(pts[:,0],pts[:,1],pts[:,2],".")
    ax.set_title(r'$X$ with $%d$ points' % num_pts, fontsize = 10)
    ax.set_zticks([1])
    ax.set_yticks([0,1])
    ax.set_xticks([0,1])
if dim_pts == 2:
    ax = fig.add_subplot(1, 2, 2)
    ax.set_title(r'$X$ with $%d$ points' % num_pts, fontsize = 10)
    plt.scatter(pts[:,0],pts[:,1], s=10)
if dim_pts == 1:
    ax = fig.add_subplot(1, 2, 2)
    ax.set_title(r'$X$ with $%d$ points' % num_pts, fontsize = 10)
    plt.scatter(pts,0*pts, s=10)


# diagrams
ax = fig.add_subplot(1, 2, 1)
ax.set_title(r'$\mathrm{Dgm}_k(X)$  $\epsilon = %.2f$' % float(ripser_threshold), fontsize = 10)
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
plt.xlabel('birth', fontsize=10)
plt.ylabel('death', fontsize=10)
adjust_spines(ax, ['left', 'bottom'])
plt.plot([0,infinity],[0,infinity], c='k', linewidth=1)
colors = ['r', 'b', 'g']
markers = ['o', 's', 'x']
labels = [r'$k = %d$' % x for x in range(len(dimensions))]

for dim in dimensions:
    #print diagrams_inf[dim]
    plt.scatter(diagrams[dim].birth, diagrams[dim].death, s=10, c=colors[dim], label=labels[dim], marker=markers[dim])
    if len(diagrams_inf[dim].index) > 0:
        plt.scatter(diagrams_inf[dim].birth, -infinity*diagrams_inf[dim].death, s=15, c=colors[dim], marker=markers[dim], label='') 
plt.legend(loc="lower right")


# barcodes
for dim in dimensions:
    print "Number of bars in dimension %d: %d" % (dim, bar_num[dim] + bar_inf_num[dim])
    fig = plt.figure()
    ax = plt.subplot("111")
    ax.set_title("%d-dimensional bars: %d finite, %d infinite" % (dim, bar_num[dim], bar_inf_num[dim]), fontsize = 10)
    # infinite bars
    if bar_inf_num[dim] > 0:
        for i in range(bar_inf_num[dim]):
                h=i+bar_num[dim]
                plot_bar([diagrams_inf[dim].birth[i],h],[-infinity*diagrams_inf[dim].death[i],h])
                plt.scatter([-infinity*diagrams_inf[dim].death[i]],[h], c='b', s=10, marker='>')
    # finite bars
    if bar_num[dim] > 0:
        for i in range(bar_num[dim]):
                plot_bar([diagrams[dim].birth[i],i],[diagrams[dim].death[i],i])
                plt.plot([diagrams[dim].death.max(),diagrams[dim].death.max()],[0,bar_num[dim]], c='r', linestyle='--', linewidth=0.5)
        #plt.xticks(list(plt.xticks()[0]) + [diagrams[dim].death.max()])

plt.show()
#'''
