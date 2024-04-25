import yt
# from yt.utilities.physical_constants import kb
# from yt import YTQuantity
# import numpy as np
# import os
import matplotlib
import matplotlib.pyplot as plt
import argparse
from mpi4py import MPI

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 12, 'font.weight' : 'bold'})

accumulated = True
USE_MPI = False # MPI functions haven't been written for this file
if USE_MPI:
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
else:
    size = int(1)
    rank = int(0)

parser = argparse.ArgumentParser(description='Calculate the radial profile of the given physical quantity.')
parser.add_argument('-n', '--number'  , nargs=1, type=int, help="number of your output file", default=[0])
parser.add_argument('-f', '--filename', nargs=1          , help='the prefix of your output file names', default=['Data_00'])
parser.add_argument('-d', '--directory'                  , help='the name of the simulation directory', default=[None])
# parser.add_argument('-d', '--directory', nargs=1, help='the name of the simulation directory', default=[None])
args = parser.parse_args()
if args.directory[0] == None:
    raise NameError('the simulation directory is not found!!')
num = args.number[0]
file_name_prefix = args.filename[0]
directory = args.directory.split()[0]
count = str(num).zfill(4)
filename = '../' + directory + '/' + file_name_prefix + count
# figure_name = 'radial_profile/density_all'
if accumulated:
    figure_name = 'radial_profile/accumulated_mass' + directory[:4] + '_' + count
else:
    pass
    # figure_name = 'radial_profile/mass_' + directory[:4] + '_' + count

ds = yt.load(filename)
sp = ds.sphere("c", (500,"kpc"))
rp = yt.create_profile(sp,'radius','mass',units={'radius':'kpc'},n_bins=256, weight_field=None, accumulation=accumulated)
time = '%.2f Myr' %ds.current_time.in_units('Myr')

fig = plt.figure()
ax  = fig.add_subplot(111)
ax.loglog(rp.x.value,rp["mass"].in_units("g").value)
# ax.plot(rp.x.value,rp["density"].in_units("g/cm**3").value)
if accumulated:
    pass
else:
    ax.set_xlim(1, 250)
    ax.set_ylim(1e-27, 1e-24)
ax.set_title('mass radial profile at '+time)
ax.set_xlabel(r'radius $(kpc)$')
ax.set_ylabel(r'mass $(g)$')
fig.savefig(figure_name)