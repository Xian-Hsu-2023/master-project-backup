import yt
from yt.utilities.physical_constants import kb
from yt import YTQuantity
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI
# mpiexec -n 45 python hardcastle14_lobe_property.py

dir_name = ["s057_t250_r3kb8_LB50_ref7_0.0", 
            "s056_t250_r3kb8_LB50_2A05_ref7_0.001",
            "s055_t250_r3kb8_LB50_2A05_ref7_0.01"]
start_ID = int(1)
max_ID = int(50)
interval = int(1)
number = int((max_ID-start_ID)/interval+1)
time = np.zeros(number)
length = np.zeros((len(dir_name), number))
volume = np.zeros((len(dir_name), number))
recv_tim = np.zeros_like(time)
recv_len = np.zeros_like(length)
recv_vol = np.zeros_like(volume)
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

for i in range(len(dir_name)):
    for j in range(number):
        if (i*len(dir_name)+j) % size == rank:
            ds = yt.load('../'+dir_name[i]+'/Data_'+str(j*interval+start_ID).zfill(6))
            sp = ds.sphere("c", (400,"kpc"))
            lobe = ds.cut_region(sp, ['obj["gas", "cooling_time"].in_units("Gyr") > 35']) # larger than 25 Gyr is necessary
            length[i, j] = lobe.argmax(("gas", "x"))[0] / YTQuantity(1, 'kpc') - 4000
            volume[i, j] = lobe.sum("cell_volume") / YTQuantity(1, 'kpc')**3
            if i == 0:
                time[j] = ds.current_time.in_units('Myr')

comm.Barrier()
comm.Reduce(length, recv_len)
comm.Reduce(volume, recv_vol)
comm.Reduce(time, recv_tim)

if rank==0:
    fig = plt.figure(figsize=(6,12))
    ax1 = fig.add_subplot(311)
    ax1.set_title('lobe property')
    # ax1.loglog(recv_tim,recv_len[0])
    for i in range(len(dir_name)):
        ax1.loglog(recv_tim,recv_len[i])
    ax1.set_xlim(1, 250)
    ax1.set_ylim(1, 400)
    ax1.set_ylabel('lobe length (kpc)')
    
    ax2 = fig.add_subplot(312)
    # ax2.loglog(recv_tim, recv_vol[0])
    for i in range(len(dir_name)):
        ax2.loglog(recv_tim, recv_vol[i])
    ax2.set_xlim(1, 250)
    ax2.set_ylabel(r'volume ($kpc^3$)')
    
    ax3 = fig.add_subplot(313)
    # ax3.loglog(recv_tim, recv_vol[0]/recv_len[0]**3)
    for i in range(len(dir_name)):
        ax3.loglog(recv_tim, recv_vol[i]/recv_len[i]**3)
    ax3.set_xlim(1, 250)
    ax3.set_xlabel('time (Myr)')
    ax3.set_ylabel(r'$volume/length^3$')
    
    for ax in [ax1, ax2, ax3]:
        ax.legend([r'$f_B=0$', r'$f_B=0.001$', r'$f_B=0.01$'])
    
    fig.savefig("hardcastle14_plots/lobe_property_s057.png")