import yt
from yt.utilities.physical_constants import kb
from yt import YTQuantity
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI
# mpiexec -n 45 python hardcastle14_plasmabeta.py

dir_name = ["s057_t250_r3kb8_LB50_ref7_0.0"]
start_ID = int(1)
max_ID = int(60)
interval = int(1)
number = int((max_ID-start_ID)/interval+1)

# time = np.zeros(number)
length     = np.zeros((len(dir_name), number))
MagOverThe = np.zeros((len(dir_name), number))
recv_len        = np.zeros_like(length)
recv_MagOverThe = np.zeros_like(MagOverThe)

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
def _thermal_energy(field, data):
    ethe = data["gas", "thermal_energy_density"]
    vol  = data["gas", "cell_volume"]
    return ethe*vol
def _magnetic_energy(field, data):
    emag = data["gas", "magnetic_energy_density"]
    vol  = data["gas", "cell_volume"]
    return emag*vol
yt.add_field(("gas", "thermal_energy") , function=_thermal_energy , sampling_type="local", units="erg")
yt.add_field(("gas", "magnetic_energy"), function=_magnetic_energy, sampling_type="local", units="erg")

for i in range(len(dir_name)):
    for j in range(number):
        if (i*len(dir_name)+j) % size == rank:
            ds = yt.load('../'+dir_name[i]+'/Data_'+str(j*interval+start_ID).zfill(6))
            sp = ds.sphere("c", (400,"kpc"))
            lobe = ds.cut_region(sp, ['obj["gas", "cooling_time"].in_units("Gyr") > 35']) # larger than 25 Gyr is necessary
            lobe_ethe = lobe.sum("thermal_energy")
            lobe_emag = lobe.sum("magnetic_energy")
            MagOverThe[i, j] = lobe_emag / lobe_ethe
            length[i, j] = lobe.argmax(("gas", "x"))[0] / YTQuantity(1, 'kpc') - 4000

comm.Barrier()
comm.Reduce(length, recv_len)
comm.Reduce(MagOverThe, recv_MagOverThe)

if rank==0:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.loglog(recv_len[0], recv_MagOverThe[0])
    ax.set_xlabel('Lobe length (kpc)')
    ax.set_ylabel('Magnetic energy / thermal energy')
    fig.savefig('hardcastle14_plots/plasmabeta_test.png')