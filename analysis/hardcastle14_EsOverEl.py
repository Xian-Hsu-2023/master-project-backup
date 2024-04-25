import yt
from yt.utilities.physical_constants import kb
from yt import YTQuantity
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI
# mpiexec -n 45 python hardcastle14_EsOverEl.py

dir_name = ["s057_t250_r3kb8_LB50_ref7_0.0"]
start_ID = int(1)
max_ID = int(60)
interval = int(1)
number = int((max_ID-start_ID)/interval+1)
length     = np.zeros((len(dir_name), number))
es_over_el = np.zeros((len(dir_name), number))
recv_len        = np.zeros_like(length)
recv_es_over_el = np.zeros_like(es_over_el)
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 12, 'font.weight' : 'bold'})

def _tCool(field, data): 
    mu = 0.61 # assuming fully ionized gas 
    mue = 1.18 
    mp = 1.67e-24*YTQuantity(1.,"g") 
    ne = data["gas", "density"]/mue/mp 
    n  = data["gas", "density"]/mu/mp 
    T  = data["gas", "temperature"].in_units('K') 
    return 1.5*n/ne**2.*kb*T**0.5/YTQuantity(3e-27, 'erg/s*cm**3/K**0.5')
yt.add_field(("gas", "cooling_time"), function = _tCool, sampling_type="local", units="Gyr") 
def _total_energy(field, data):
    etot = data["gas", "total_energy_density"]
    vol  = data["gas", "cell_volume"]
    return etot*vol
yt.add_field(("gas", "total_energy")   , function=_total_energy, sampling_type="local", units="erg")
def _hardcastle14_radial_mach_number(field, data):
    rx = data["gas", "x"] - YTQuantity(4000, 'kpc')
    ry = data["gas", "y"] - YTQuantity(4000, 'kpc')
    rz = data["gas", "z"] - YTQuantity(4000, 'kpc')
    vx = data["gas", "velocity_x"]
    vy = data["gas", "velocity_y"]
    vz = data["gas", "velocity_z"]
    cs = data["gas", "sound_speed"]
    result = (vx*rx+vy*ry+vz*rz)/(rx**2+ry**2+rz**2)**0.5/cs
    return abs(result)
yt.add_field(("gas", "hardcastle14_radial_mach_number"), function=_hardcastle14_radial_mach_number, sampling_type="local", units="")

for i in range(len(dir_name)):
    for j in range(number):
        if (i*len(dir_name)+j) % size == rank:
            ds = yt.load('../'+dir_name[i]+'/Data_'+str(j*interval+start_ID).zfill(6))
            sp = ds.sphere("c", (400,"kpc"))
            lobe    = ds.cut_region(sp, ['obj["gas", "cooling_time"].in_units("Gyr") > 35']) # larger than 25 Gyr is necessary
            notlobe = ds.cut_region(sp, ['obj["gas", "cooling_time"].in_units("Gyr") < 35']) # larger than 25 Gyr is necessary
            shock   = ds.cut_region(notlobe, ['obj["gas", "hardcastle14_radial_mach_number"] > 0.1'])
            es = shock.sum("total_energy")
            el =  lobe.sum("total_energy")
            es_over_el[i, j] = es/el
            length[i, j] = lobe.argmax(("gas", "x"))[0] / YTQuantity(1, 'kpc') - 4000
            # if i == 0:
            #     time[j] = ds.current_time.in_units('Myr')

comm.Barrier()
comm.Reduce(length, recv_len)
comm.Reduce(es_over_el, recv_es_over_el)

if rank==0:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(recv_len[0], recv_es_over_el[0], label='')
    ax.set_xlabel('Lobe length (kpc)')
    ax.set_ylabel('Shocked region energy/Lobe energy')
    fig.savefig('hardcastle14_plots/EsOverEl_test.png')