import yt
# from yt.utilities.physical_constants import kb
# from yt import YTQuantity
import numpy as np
# import os
import matplotlib
import matplotlib.pyplot as plt
import argparse
from mpi4py import MPI

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 12, 'font.weight' : 'bold'})

accumulated = False
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
    figure_name = 'radial_profile/accumulate_density_' + directory[:4] + '_' + count
else:
    figure_name = 'radial_profile/density_' + directory[:4] + '_' + count

ds = yt.load(filename)
sp = ds.sphere("c", (500,"kpc"))
# rp = yt.create_profile(sp,'radius','density',units={'radius':'kpc'},n_bins=256, weight_field=None, accumulation=True)
rp = yt.create_profile(sp,'radius','density',units={'radius':'kpc'},n_bins=256, weight_field=("gas", "cell_volume"), accumulation=False)
# rp = yt.create_profile(sp,'radius','density',units={'radius':'kpc'},n_bins=256, weight_field=("gas", "cell_volume"), accumulation=accumulated)
time = '%.2f Myr' %ds.current_time.in_units('Myr')

fig = plt.figure()
ax  = fig.add_subplot(111)
ax.loglog(rp.x.value,rp["density"].in_units("g/cm**3").value)

filename = '/lfs/data/xianhsu/flash/mhd_injectiontest131_toroidal_kincoef175_ONeill_namb_inj10_ref8_0.001/crbub_hdf5_plt_cnt_0000'
ds = yt.load(filename)
sp = ds.sphere("c", (500,"kpc"))
rp = yt.create_profile(sp,'radius','density',units={'radius':'kpc'},n_bins=256, weight_field=("gas", "cell_volume"), accumulation=False)
ax.loglog(rp.x.value,rp["density"].in_units("g/cm**3").value)

Msun = 1.9891e33 # g
Mpc = 3.0857e24 # cm
Kpc = 3.0857e21 # cm
M = 8.5e14*Msun
Rvir = 2.381*Mpc
conc = 6.81
A_NFW = np.log(1+conc) - conc/(1+conc)
rho_halo = M/(4/3*np.pi*Rvir**3)
r = np.logspace(-3, 0, num=200) * 0.5 * Mpc
x = r/Rvir
hubble = 0.7
m = M*1.0e-15/hubble
fgas = ((0.72/hubble)**1.5)*(0.125+0.037*np.log10(m))
rho = rho_halo / (3 * A_NFW * x * (1/conc+x)**2)
ax.loglog(r/Kpc, rho, '-')

if accumulated:
    pass
else:
    ax.set_xlim(1, 500)
    # ax.set_ylim(1e-27, 1e-24)
ax.set_title('density radial profile at '+time)
ax.set_xlabel(r'radius $(kpc)$')
ax.set_ylabel(r'density $(g\, cm^{-3})$')
ax.legend(['gamer', 'flash', 'analytical'])
fig.savefig(figure_name)