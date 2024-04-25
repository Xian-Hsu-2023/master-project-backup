import yt
from yt.utilities.physical_constants import kb
from yt import YTQuantity
# import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.colors import LogNorm
import os
import argparse
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

parser = argparse.ArgumentParser(description='Plot the cooling time with given plot file numbers, file name, region, and direction.')
parser.add_argument('-n', '--number' , nargs=1, type=int, help="number of the output files", default=[1])
parser.add_argument('-f', '--filename', nargs=1, help='the name of your output files', default=['Data_'])
parser.add_argument('-r', '--region', nargs=1, type=int, help='the size (kpc) for your image', default=[100])
parser.add_argument('-d', '--direction', nargs=1, help='the direction (x, y, or z) of the plotted image', default=['y'])
parser.add_argument('-v', '--variable', nargs=1, help='the variable that you want to plot', default=['density'])
# parser.add_argument('-dir', '--directory', nargs=1, help='the directory of the output files', default='')
args = parser.parse_args()
num = args.number[0]
file_name_type = args.filename[0]
region = args.region[0]
direction = args.direction[0]
variable = args.variable[0]
# directory = args.directory
titles = [r"$ref7, R_j=2\, kpc, A=0.5, f_B=0.1$", r"$f_B=0.01$"]
# titles = [r"L_B_{MAX}=50kpc", r"L_B_{MAX}=100kpc", r"L_B_{MAX}=500kpc"]
# file_directory = ['s042_regular_r3000b8_ambLB50_ref7_0.1', 's041_regular_r3000b8_ambLB100_ref7_0.1', 's039_regular_r3000b8_amb_ref7_0.1']
file_directory = ['s054_t250_r3kb8_LB50_2A05_ref7_0.1', 's055_t250_r3kb8_LB50_2A05_ref7_0.01']
directory = 'figure/magtest' + file_directory[0][1:4] + '_' + variable + '_' + direction + '_range' + str(region)
command = 'mkdir ' + directory
if rank==0:
    os.system(command)
# else:
#     print('Hello from rank %.f2' %rank, flush=True)
    
count = [None]*1000
for i in range(1000):
    count[i] = str(i).zfill(6)

def _tCool(field, data): 
    mu = 0.61 # assuming fully ionized gas 
    mue = 1.18 
    mp = 1.67e-24*YTQuantity(1.,"g") 
    ne = data["gas", "density"]/mue/mp 
    n  = data["gas", "density"]/mu/mp 
    T  = data["gas", "temperature"].in_units('K') 
    return 1.5*n/ne**2.*kb*T**0.5/YTQuantity(3e-27, 'erg/s*cm**3/K**0.5')
yt.add_field(("gas", "cooling_time"), function = _tCool, sampling_type="local", units="Gyr")
def _faraday_rotation_x(field, data):
    coef = 0.812*3.2407792896664E-19 
    mue = 1.18 
    mp = 1.67e-24*YTQuantity(1.,"g") 
    ne = data["density"]/mue/mp
    return coef * ne * data["magnetic_field_x"] * YTQuantity(1e6, 'rad/m**2*cm**3/gauss/cm')
yt.add_field(("gas", "faraday_rotation_x"), function = _faraday_rotation_x, sampling_type="local", units="rad/m**2/cm")
def _faraday_rotation_y(field, data):
    coef = 0.812*3.2407792896664E-19 
    mue = 1.18 
    mp = 1.67e-24*YTQuantity(1.,"g") 
    ne = data["density"]/mue/mp
    return coef * ne * data["magnetic_field_y"] * YTQuantity(1e6, 'rad/m**2*cm**3/gauss/cm')
yt.add_field(("gas", "faraday_rotation_y"), function = _faraday_rotation_y, sampling_type="local", units="rad/m**2/cm")
def _faraday_rotation_z(field, data):
    coef = 0.81*3.2407792896664E-19 
    mue = 1.18 
    mp = 1.67e-24*YTQuantity(1.,"g") 
    ne = data["density"]/mue/mp
    return coef * ne * data["magnetic_field_z"] * YTQuantity(1e6, 'rad/m**2*cm**3/gauss/cm')
yt.add_field(("gas", "faraday_rotation_z"), function = _faraday_rotation_z, sampling_type="local", units="rad/m**2/cm")


for i in range(num):
    if (i % size == rank):
        # problem: can't control figure size with this parameter (solved)
        figx = 6 * len(file_directory)
        plt.rcParams.update({'figure.figsize': [figx, 6.]})
        fig = plt.figure()
        plt.rcParams.update({'font.size': 12})
        grid = ImageGrid(
            fig,
            (0.1, 0.1, 0.8, 0.85),
            # 111,
            nrows_ncols=(1, len(file_directory)),
            axes_pad=0.01,
            label_mode="L",
            aspect=True,
            share_all=False,
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
            if variable=='faraday_rotation_x' or variable=='faraday_rotation_y' or variable=='faraday_rotation_z':
                p = yt.ProjectionPlot(ds, direction, variable, width=(region, 'kpc'))
                p.set_zlim(variable, -1e5, 1e5)
                p.set_log(variable, True, linthresh=1e3)
            elif variable=='particle_mass':            
                p = yt.ParticlePlot(ds, ("io","particle_position_x"), ("io","particle_position_y"), z_fields=("io","particle_mass"), width=(region, "kpc"))
            else:
                p = yt.SlicePlot(ds, direction, variable, width=(region, "kpc"))
            # p.annotate_grids()
            if variable=='density':
                p.set_zlim(("gas", variable), 1e-27,2e-25)
            elif variable=='kinetic_energy_density':
                p.set_zlim(("gas", variable), 1e-11,1e-16)
            elif variable=='cooling_time':
                p.set_zlim(("gas", variable), 1e-1,30)
            elif variable=='magnetic_energy_density':
                p.set_zlim(("gas", variable), 1e-15,1e-10)
            elif variable=='four_velocity_magnitude':
                p.set_zlim(("gas", variable), 1e4,1e7)
            elif variable=='magnetic_field_magnitude':
                p.set_zlim(("gas", variable), 1e-6, 5e-5)
                p.annotate_magnetic_field(headlength=3)
            elif variable=='temperature':
                p.set_zlim(('gas', variable), 1e7, 1e12)
            elif variable=='ParDens':
                p.set_zlim((variable), 1e-6, 5e2)
            # p.set_zlim(("gas", "density"), 1e-25,2e-24)
            # p.annotate_grids()
            
            if variable=='particle_mass':
                plot = p.plots[('io', variable)]
            else:
                plot = p.plots[(variable)]
            # plot = p.plots[("gas", variable)]
            # plot.figure = fig
            plot.axes = grid[j].axes
            plot.cax = grid.cbar_axes[j]
            p.render()
            grid[j].axes.set_title(titles[j])
            # plots[j].set_cmap("gnuplot")
            # plots[j].set_clim((1,2e2))
            # plots[j].annotate_grids()
            
        textx = 0.2
        # textx = 0.15+0.05*len(file_directory)
        texty = 0.925
        grid[0].axes.text(textx, texty, 'time = %.2f Myr' %ds.current_time.in_units('Myr'), horizontalalignment='center', 
                          color="black", bbox=dict(facecolor='white', alpha=0.5), verticalalignment='center', transform=grid[0].axes.transAxes)
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