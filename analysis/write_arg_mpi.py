#     usage: python write_arg_mpi.py -n $file_number -f Data_ -d $directory_name
#       MPI: mpiexec -n 10 python write_arg_mpi.py -n $file_number -f Data_ -d $directory_name
# first use: execute "mkdir text_file"
import yt
from yt.units import Msun, Mpc, dyne, cm, g
from yt.utilities.physical_constants import kb
from yt import YTQuantity
import numpy as np
import os
import argparse
from mpi4py import MPI
# MPI switch
USE_MPI = True
if USE_MPI:
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
else:
    size = int(1)
    rank = int(0)

parser = argparse.ArgumentParser(description='Sum the energies from the simulation box, and wrtie them into a text file.')
parser.add_argument('-n', '--number' , nargs=1, type=int, help="number of your output files", default=[101])
parser.add_argument('-f', '--filename', nargs=1, help='the prefix of your output file names', default=['Data_'])
parser.add_argument('-d', '--directory', nargs=1, help='the name of the simulation directory', default=[None])
args = parser.parse_args()
if args.directory[0] == None:
    raise NameError('the simulation directory is not found!!')
num = args.number[0]
file_name_prefix = args.filename[0]
directory = args.directory[0]
# mkdir text_file
textfile_name = 'text_file/Record__' + directory + '.txt'
count = [None]*num
for i in range(num):
    count[i] = str(i).zfill(6)
G = 6.67430e-8*dyne*cm**2/g**2
M = 8.5e14*Msun # Perseus cluster
r_vir = 2.44*Mpc
# M = 1e14*Msun # ZuHone model
# r_vir = 1.201*Mpc
c = 4
# c = 6.81
t        = np.zeros(num)
E_k_dif  = np.zeros(num)
E_th_dif = np.zeros(num)
E_g_dif  = np.zeros(num)
E_m_dif  = np.zeros(num)
E_dif    = np.zeros(num)
if USE_MPI:
    t_sum        = np.zeros(num)
    E_k_dif_sum  = np.zeros(num)
    E_th_dif_sum = np.zeros(num)
    E_g_dif_sum  = np.zeros(num)
    E_m_dif_sum  = np.zeros(num)
    E_dif_sum    = np.zeros(num)

# add the user-defined fields (since latest yt version may not support "from yt.mods import *" method)
# (if encounter errors with default fields in data, check the field list in yt for your own simulation)
def _tCool(field, data):   # cooling time of the ICM
    mu = 0.61              # assuming fully ionized gas 
    mue = 1.18 
    mp = 1.67e-24*YTQuantity(1.,"g") 
    ne = data["gas", "density"]/mue/mp 
    n  = data["gas", "density"]/mu/mp 
    T  = data["gas", "temperature"].in_units('K') 
    return 1.5*n/ne**2.*kb*T**0.5/YTQuantity(3e-27, 'erg/s*cm**3/K**0.5')
def _E_k(field, data):     # kinetic energy
    return data["density"]/2*data["velocity_magnitude"]**2*data["cell_volume"] 
def _E_th(field, data):    # thermal energy
    return data["thermal_energy_density"]*data["cell_volume"] 
# def _phi_NFW(field, data): # gravitational potential
#     return -G*M/data["radius"]*np.log(1+data["radius"]/(r_vir/c))/(np.log(1+c)-c/(1+c)) 
# def _E_g(field, data):     # gravitational energy
#     return data["density"]*data["phi_NFW"]*data["cell_volume"] 
def _E_g(field, data):     # gravitational energy
    return (-data["mass"])*data["gravitational_potential"]
def _E_m(field, data):     # magnetic energy
    return data["magnetic_energy_density"]*data["cell_volume"] 
yt.add_field(("gas", "cooling_time"), function = _tCool, sampling_type="local", units="Gyr")
yt.add_field(("gas", "E_k"), function = _E_k, sampling_type="local", units="erg")
yt.add_field(("gas", "E_th"), function = _E_th, sampling_type="local", units="erg")
# yt.add_field(("gas", "phi_NFW"), function = _phi_NFW, sampling_type="local", units="erg/g")
yt.add_field(("gas", "E_g"), function = _E_g, sampling_type="local", units="erg")
yt.add_field(("gas", "E_m"), function = _E_m, sampling_type="local", units="erg")

file_name = '../' + directory + '/' + file_name_prefix + count[0]
ds_zero = yt.load(file_name)
ad = ds_zero.all_data()
E_k_initial  = ad["gas", "E_k"].sum()
E_th_initial = ad["gas", "E_th"].sum()
E_g_initial  = ad["gas", "E_g"].sum()
E_m_initial  = ad["gas", "E_m"].sum()

for i in range(num):
    if (i % size) == rank:
        file_name = '../' + directory + '/' + file_name_prefix + count[i]
        ds = yt.load(file_name)
        ad = ds.all_data()
        t[i] = ds.current_time.in_units('Myr')
        E_k_dif[i]  = ad["gas", "E_k"].sum()  - E_k_initial
        E_th_dif[i] = ad["gas", "E_th"].sum() - E_th_initial
        E_g_dif[i]  = ad["gas", "E_g"].sum()  - E_g_initial
        E_m_dif[i]  = ad["gas", "E_m"].sum()  - E_m_initial
        E_dif[i]    = E_k_dif[i] + E_th_dif[i] + E_g_dif[i] + E_m_dif[i]

if USE_MPI:
    comm.Reduce(t       , t_sum       , op = MPI.SUM, root = 0)
    comm.Reduce(E_k_dif , E_k_dif_sum , op = MPI.SUM, root = 0)
    comm.Reduce(E_th_dif, E_th_dif_sum, op = MPI.SUM, root = 0)
    comm.Reduce(E_g_dif , E_g_dif_sum , op = MPI.SUM, root = 0)
    comm.Reduce(E_m_dif , E_m_dif_sum , op = MPI.SUM, root = 0)
    comm.Reduce(E_dif   , E_dif_sum   , op = MPI.SUM, root = 0)
    t        = np.copy(t_sum)
    E_k_dif  = np.copy(E_k_dif_sum)
    E_th_dif = np.copy(E_th_dif_sum)
    E_g_dif  = np.copy(E_g_dif_sum)
    E_m_dif  = np.copy(E_m_dif_sum)
    E_dif    = np.copy(E_dif_sum)

if rank==0:
    f = open(textfile_name, "w")
    f.write(f"{'time(Myr)'.ljust(13)} {'kinetic(erg)'.ljust(13)} {'thermal'.ljust(13)} {'gravitational'.ljust(13)} {'magnetic'.ljust(13)} {'total'.ljust(13)} \n")
    # f.write(f"{'time(Myr)':.13s} {'kinetic(erg)':.13s} {'thermal':.13s} {'gravitational':.13s} {'magnetic':.13s} {'total':.13s} \n")
    for i in range(num):
        f.write(f"{t[i]:.7e} {E_k_dif[i]:.7e} {E_th_dif[i]:.7e} {E_g_dif[i]:+.6e} {E_m_dif[i]:.7e} {E_dif[i]:.7e} \n")
        # f.write("%.3e %.3e %.3e %.3e %.3e %.3e \n" %(t[i], E_k_dif[i], E_th_dif[i], E_g_dif[i], E_m_dif[i], E_dif[i]))
    f.close()