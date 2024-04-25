import yt
from yt.utilities.physical_constants import kb, pi
# from yt.visualization.base_plot_types import get_multi_plot
from yt import YTQuantity
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.colors import LogNorm
import os
import argparse

parser = argparse.ArgumentParser(description='Plot the cooling time with given plot file numbers, file name, region, and direction.')
parser.add_argument('-n', '--number' , nargs=1, type=int, help="number of the output files", default=[100])
parser.add_argument('-f', '--filename', nargs=1, help='the name of your output files', default=['Data_'])
parser.add_argument('-r', '--region', nargs=1, type=int, help='the size (kpc) for your image', default=[100])
parser.add_argument('-d', '--direction', nargs=1, help='the direction (x, y, or z) of the plotted image', default=['x'])
# parser.add_argument('-dir', '--directory', nargs=1, help='the directory of the output files', default='')
args = parser.parse_args()
num = args.number[0]
file_name_type = args.filename[0]
region = args.region[0]
direction = args.direction[0]
# directory = args.directory
titles=['lv2','lv3','lv4']
# titles=['lv1','lv2','lv3']
file_directory = ['singletest003_ref2',
                  'singletest002_ref3',
                  'singletest001_ref4']
# file_directory = ['singletest_ref1',
#                   'singletest_ref2',
#                   'singletest_ref3']
count = [None]*1000
for i in range(1000):
    count[i] = str(i).zfill(6)
# plots = [None]*3
# orient = 'horizontal'

def _tCool(field, data): 
    mu = 0.61 # assuming fully ionized gas 
    mue = 1.18 
    mp = 1.67e-24*YTQuantity(1.,"g") 
    ne = data["gas", "density"]/mue/mp 
    n  = data["gas", "density"]/mu/mp 
    T  = data["gas", "temperature"].in_units('K') 
    return 1.5*n/ne**2.*kb*T**0.5/YTQuantity(3e-27, 'erg/s*cm**3/K**0.5')
yt.add_field(("gas", "cooling_time"), function = _tCool, sampling_type="local", units="Gyr")

directory = 'coolingtime_'+direction+'_range'+str(region)
command = 'mkdir '+directory
os.system(command)
for i in range(num):    
    fig = plt.figure()
    plt.rcParams.update({'font.size': 10})
    grid = AxesGrid(
        fig,
        (0.1, 0.1, 0.8, 0.9),
        nrows_ncols=(1, 3),
        axes_pad=0.01,
        label_mode="L",
        share_all=True,
        cbar_location="right",
        cbar_mode="single",
        cbar_size="3%",
        cbar_pad="0%",
    )
    # titles=['lv1','lv2','lv3']
    for j in range(len(file_directory)):
        file_name = '../'+file_directory[j]+'/'+file_name_type+count[i]
        #ds = yt.load('../crbub_hdf5_plt_cnt_'+count[i])
        ds = yt.load(file_name)
        p = yt.SlicePlot(ds, direction, 'cooling_time', width=(region, "kpc"))
        p.annotate_grids()
        p.set_zlim(("gas", "cooling_time"), 1, 2e2)
        p.annotate_grids()
        plot = p.plots[("gas", "cooling_time")]
        plot.figure = fig
        plot.axes = grid[j].axes
        plot.cax = grid.cbar_axes[j]
        p.render()
        grid[j].axes.set_title(titles[j])
        # plots[j].set_cmap("gnuplot")
        # plots[j].set_clim((1,2e2))
        # plots[j].annotate_grids()
    grid[0].axes.text(0.3, 0.9, 'time = %.2f Myr' %ds.current_time.in_units('Myr'), horizontalalignment='center',verticalalignment='center', transform=grid[0].axes.transAxes)
    # axes[0][0].text(20,30,'time = %.2f Myr' %ds.current_time.in_units('Myr'),color='w',size=15)
    #axes[0][0].text(5,760,'density',color='w',size=15)
    #axes[0][1].text(5,760,'pressure',color='w',size=15)
    #axes[0][2].text(5,760,'temperature',color='w',size=15)
    try:
        fig.savefig(directory+'/'+count[i]+'.png')
    except:
        print("error happens when saving file... QQ")
    plt.close(fig)
    del fig, grid
    # del fig, axes, colorbars