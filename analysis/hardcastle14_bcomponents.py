import yt
from yt.utilities.physical_constants import kb
from yt import YTQuantity
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI
# mpiexec -n 45 python hardcastle14_bcomponents.py

dir_name = ["s057_t250_r3kb8_LB50_ref7_0.0"]
start_ID = int(1)
max_ID = int(90)
interval = int(1)
number = int((max_ID-start_ID)/interval+1)

time = np.zeros(number)
length = np.zeros((len(dir_name), number))
frac_tor = np.zeros((len(dir_name), number))
frac_lon = np.zeros((len(dir_name), number))
frac_rad = np.zeros((len(dir_name), number))

recv_tim = np.zeros_like(time)
recv_len = np.zeros_like(length)
recv_frac_tor = np.zeros_like(frac_tor)
recv_frac_lon = np.zeros_like(frac_lon)
recv_frac_rad = np.zeros_like(frac_rad)
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
def _magnetic_energy(field, data):
    eb = data["gas", "magnetic_energy_density"]
    vol= data["gas", "cell_volume"]
    return eb*vol
yt.add_field(("gas", "magnetic_energy"), function=_magnetic_energy, sampling_type="local", units="erg")
def _magnetic_energy_longitudinal(field, data):
    bx = data["gas", "magnetic_field_x"]
    vol= data["gas", "cell_volume"]
    return bx**2/(8*np.pi) * vol
def _magnetic_energy_toroidal(field, data):
    by = data["gas", "magnetic_field_y"]
    bz = data["gas", "magnetic_field_z"]
    ry = data["gas", "y"] - YTQuantity(4000, 'kpc')
    rz = data["gas", "z"] - YTQuantity(4000, 'kpc')
    b  = (by*rz-bz*ry)/(ry**2+rz**2)**0.5
    vol= data["gas", "cell_volume"]
    return b**2/(8*np.pi) * vol
def _magnetic_energy_radial(field, data):
    by = data["gas", "magnetic_field_y"]
    bz = data["gas", "magnetic_field_z"]
    ry = data["gas", "y"] - YTQuantity(4000, 'kpc')
    rz = data["gas", "z"] - YTQuantity(4000, 'kpc')
    b  = (by*ry+bz*rz)/(ry**2+rz**2)**0.5
    vol= data["gas", "cell_volume"]
    return b**2/(8*np.pi) * vol
yt.add_field(("gas", "magnetic_energy_longitudinal"), function=_magnetic_energy_longitudinal, sampling_type="local", units="erg")
yt.add_field(("gas", "magnetic_energy_toroidal")    , function=_magnetic_energy_toroidal    , sampling_type="local", units="erg")
yt.add_field(("gas", "magnetic_energy_radial")      , function=_magnetic_energy_radial      , sampling_type="local", units="erg")

for i in range(len(dir_name)):
    for j in range(number):
        if (i*len(dir_name)+j) % size == rank:
            ds = yt.load('../'+dir_name[i]+'/Data_'+str(j*interval+start_ID).zfill(6))
            sp = ds.sphere("c", (400,"kpc"))
            lobe = ds.cut_region(sp, ['obj["gas", "cooling_time"].in_units("Gyr") > 35']) # larger than 25 Gyr is necessary
            eba= lobe.sum("magnetic_energy")
            # print("finish step 1", flush=True)
            ebt= lobe.sum("magnetic_energy_toroidal")
            # print("finish step 2", flush=True)
            ebl= lobe.sum("magnetic_energy_longitudinal")
            # print("finish step 3", flush=True)
            ebr= lobe.sum("magnetic_energy_radial")
            # print(eba, ebt, ebl, ebr, (ebt+ebl+ebr)/eba)
            length[i, j] = lobe.argmax(("gas", "x"))[0] / YTQuantity(1, 'kpc') - 4000
            frac_tor[i, j] = ebt / eba
            frac_lon[i, j] = ebl / eba
            frac_rad[i, j] = ebr / eba
            if i == 0:
                time[j] = ds.current_time.in_units('Myr')

comm.Barrier()
comm.Reduce(length, recv_len)
comm.Reduce(frac_tor, recv_frac_tor)
comm.Reduce(frac_lon, recv_frac_lon)
comm.Reduce(frac_rad, recv_frac_rad)
comm.Reduce(time, recv_tim)

if rank==0:
    print(recv_tim)
    print(recv_len)
    print(recv_frac_tor)
    print(recv_frac_lon)
    print(recv_frac_rad)
    fig = plt.figure()
    ax0 = fig.add_subplot(121)
    # for i in range(len(dir_name)):
    #     ax.plot(recv_tim,recv_len[i], alpha=0.5)
    ax0.plot(recv_tim,recv_len[0])
    ax0.set_xlim(1, 250)
    ax0.set_ylim(1, 400)
    ax0.set_xlabel('time (Myr)')
    ax0.set_ylabel('jet length (kpc)')
    ax1 = fig.add_subplot(122)
    ax1.plot(recv_len[0], recv_frac_tor[0], label='toroidal')
    ax1.plot(recv_len[0], recv_frac_lon[0], label='longitudinal')
    ax1.plot(recv_len[0], recv_frac_rad[0], label='radial')
    ax1.set_xlabel(r'jet length (kpc)')
    ax1.set_ylabel(r'$\frac{E_{B_{comp}}}{E_B}$')
    ax1.legend()
    fig.savefig("hardcastle14_plots/bcomp_test.png")

