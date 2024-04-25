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
directory = args.directory[0]
count = str(num).zfill(4)
filename = '/lfs/data/xianhsu/flash/mhd_injectiontest131_toroidal_kincoef175_ONeill_namb_inj10_ref8_0.001/crbub_hdf5_plt_cnt_0000'
# filename = '../' + directory + '/' + file_name_prefix + count
if accumulated:
    figure_name = 'radial_profile/accumulated_mass_flash'
else:
    pass
    # figure_name = 'radial_profile/density_flash'
# figure_name = 'radial_profile/density_' + directory[:4] + '_' + count

ds = yt.load(filename)
sp = ds.sphere("c", (500,"kpc"))
# rp = yt.create_profile(sp,'radius','density',units={'radius':'kpc'},n_bins=256, weight_field=None, accumulation=accumulated)
rp = yt.create_profile(sp,'radius','mass',units={'radius':'kpc'},n_bins=256, weight_field=None, accumulation=accumulated)
time = '%.2f Myr' %ds.current_time.in_units('Myr')

fig = plt.figure()
ax  = fig.add_subplot(111)
ax.loglog(rp.x.value,rp["mass"].in_units("g").value)

filename = '/lfs/data/xianhsu/gamer/bin/mhd_injection/s005_CompareToFLASH_ref4_0.001/Data_000000'
ds = yt.load(filename)
sp = ds.sphere("c", (500,"kpc"))
rp = yt.create_profile(sp,'radius','mass',units={'radius':'kpc'},n_bins=256, weight_field=None, accumulation=accumulated)
ax.loglog(rp.x.value,rp["mass"].in_units("g").value)
filename = '/lfs/data/xianhsu/gamer/bin/mhd_injection/s008_ref4_0.001/Data_000000'
ds = yt.load(filename)
sp = ds.sphere("c", (500,"kpc"))
rp = yt.create_profile(sp,'radius','mass',units={'radius':'kpc'},n_bins=256, weight_field=None, accumulation=accumulated)
ax.loglog(rp.x.value,rp["mass"].in_units("g").value)
ax.legend(['flash', 'gamer, c=4.0', 'gamer, c=1'])

# ax.plot(rp.x.value,rp["density"].in_units("g/cm**3").value)
if accumulated:
    pass
else:
    ax.set_xlim(1, 500)
    ax.set_ylim(1e-27, 1e-24)
ax.set_title('mass radial profile at '+time)
ax.set_xlabel(r'radius $(kpc)$')
ax.set_ylabel(r'mass $(g)$')
fig.savefig(figure_name)