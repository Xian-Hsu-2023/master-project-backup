import yt
from yt.utilities.physical_constants import kb
from yt import YTQuantity
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI
# mpiexec -n 45 python hardcastle14_energy.py

dir_name = ["s057_t250_r3kb8_LB50_ref7_0.0"]
start_ID = int(1)
max_ID = int(60)
interval = int(1)
number = int((max_ID-start_ID)/interval+1)
time = np.zeros(number)
etot = np.zeros((len(dir_name), number))
lobe_etot = np.zeros((len(dir_name), number))
lobe_ethe  = np.zeros((len(dir_name), number))
lobe_ekin  = np.zeros((len(dir_name), number))
lobe_emag  = np.zeros((len(dir_name), number))
shock_etot = np.zeros((len(dir_name), number))
shock_ethe  = np.zeros((len(dir_name), number))
shock_ekin  = np.zeros((len(dir_name), number))
shock_emag  = np.zeros((len(dir_name), number))
lobe_emag_tor = np.zeros((len(dir_name), number))
lobe_emag_lon = np.zeros((len(dir_name), number))
lobe_emag_rad = np.zeros((len(dir_name), number))

recv_etot      = np.zeros_like(etot)
recv_time       = np.zeros_like(time)
recv_lobe_etot  = np.zeros_like(lobe_etot)
recv_lobe_ethe  = np.zeros_like(lobe_etot)
recv_lobe_ekin  = np.zeros_like(lobe_etot)
recv_lobe_emag  = np.zeros_like(lobe_etot)
recv_shock_etot = np.zeros_like(lobe_etot)
recv_shock_ethe = np.zeros_like(lobe_etot)
recv_shock_ekin = np.zeros_like(lobe_etot)
recv_shock_emag = np.zeros_like(lobe_etot)
recv_lobe_emag_tor  = np.zeros_like(lobe_etot)
recv_lobe_emag_lon  = np.zeros_like(lobe_etot)
recv_lobe_emag_rad  = np.zeros_like(lobe_etot)

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
def _total_energy(field, data):
    etot = data["gas", "total_energy_density"]
    vol  = data["gas", "cell_volume"]
    return etot*vol
def _thermal_energy(field, data):
    ethe = data["gas", "thermal_energy_density"]
    vol  = data["gas", "cell_volume"]
    return ethe*vol
def _kinetic_energy(field, data):
    ekin = data["gas", "kinetic_energy_density"]
    vol  = data["gas", "cell_volume"]
    return ekin*vol
def _magnetic_energy(field, data):
    emag = data["gas", "magnetic_energy_density"]
    vol  = data["gas", "cell_volume"]
    return emag*vol
yt.add_field(("gas", "total_energy")   , function=_total_energy   , sampling_type="local", units="erg")
yt.add_field(("gas", "thermal_energy") , function=_thermal_energy , sampling_type="local", units="erg")
yt.add_field(("gas", "kinetic_energy") , function=_kinetic_energy , sampling_type="local", units="erg")
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
            lobe    = ds.cut_region(sp, ['obj["gas", "cooling_time"].in_units("Gyr") > 35']) # larger than 25 Gyr is necessary
            notlobe = ds.cut_region(sp, ['obj["gas", "cooling_time"].in_units("Gyr") < 35']) # larger than 25 Gyr is necessary
            shock   = ds.cut_region(notlobe, ['obj["gas", "hardcastle14_radial_mach_number"] > 0.1'])
            lobe_etot[i, j] = lobe.sum("total_energy")
            lobe_ethe[i, j] = lobe.sum("thermal_energy")
            lobe_ekin[i, j] = lobe.sum("kinetic_energy")
            # lobe_emag[i, j] = lobe.sum("magnetic_energy")
            shock_etot[i, j] = shock.sum("total_energy")
            shock_ethe[i, j] = shock.sum("thermal_energy")
            shock_ekin[i, j] = shock.sum("kinetic_energy")
            shock_emag[i, j] = shock.sum("magnetic_energy")
            lobe_emag_tor[i, j] = lobe.sum("magnetic_energy_toroidal")
            lobe_emag_lon[i, j] = lobe.sum("magnetic_energy_longitudinal")
            lobe_emag_rad[i, j] = lobe.sum("magnetic_energy_radial")
            lobe_emag[i, j] = lobe_emag_tor[i, j] + lobe_emag_lon[i, j] + lobe_emag_rad[i, j]
            etot[i, j] = lobe_etot[i, j] + shock_etot[i, j]
            if i == 0:
                time[j] = ds.current_time.in_units('Myr')

comm.Barrier()
comm.Reduce(etot, recv_etot)
comm.Reduce(time, recv_time)
comm.Reduce(lobe_etot, recv_lobe_etot)
comm.Reduce(lobe_ethe, recv_lobe_ethe)
comm.Reduce(lobe_ekin, recv_lobe_ekin)
comm.Reduce(lobe_emag, recv_lobe_emag)
comm.Reduce(shock_etot, recv_shock_etot)
comm.Reduce(shock_ethe, recv_shock_ethe)
comm.Reduce(shock_ekin, recv_shock_ekin)
comm.Reduce(shock_emag, recv_shock_emag)
comm.Reduce(lobe_emag_tor, recv_lobe_emag_tor)
comm.Reduce(lobe_emag_lon, recv_lobe_emag_lon)
comm.Reduce(lobe_emag_rad, recv_lobe_emag_rad)

if rank==0:
    MYR = 86400 * 365.25 * 1e6
    eexp = np.array([5e45*recv_time[i]*MYR for i in range(len(recv_time))])
    fig = plt.figure(figsize=(12,4))
    
    ax1 = fig.add_subplot(121)
    ax1.plot(recv_time, eexp              , c='cyan'          , label='expected total energy')
    ax1.plot(recv_time, recv_etot[0]      , c='blue' , lw=5   , label='total (lobe+shock)')
    ax1.plot(recv_time, recv_lobe_etot[0] , c='red'  , lw=2   , label='lobe total')
    ax1.plot(recv_time, recv_shock_etot[0], c='green', lw=2   , label='shock total')
    ax1.plot(recv_time, recv_lobe_ethe[0] , c='red'  , ls='-' , label='lobe thermal')
    ax1.plot(recv_time, recv_lobe_ekin[0] , c='red'  , ls='--', label='lobe kinetic')
    ax1.plot(recv_time, recv_lobe_emag[0] , c='red'  , ls='-.', label='lobe magnetic')
    ax1.plot(recv_time, recv_shock_ethe[0], c='green', ls='-' , label='shock thermal')
    ax1.plot(recv_time, recv_shock_ekin[0], c='green', ls='--', label='shock kinetic')
    ax1.plot(recv_time, recv_shock_emag[0], c='green', ls='-.', label='shock magnetic')
    ax1.legend(framealpha=0.4)
    ax1.set_xlabel('Time (Myr)')
    ax1.set_ylabel('Energy (erg)')
    
    ax2 = fig.add_subplot(122)
    ax2.plot(recv_time,recv_lobe_emag[0]    , c='blue'   , label='total')
    ax2.plot(recv_time,recv_lobe_emag_tor[0], c='red'    , label='toroidal')
    ax2.plot(recv_time,recv_lobe_emag_lon[0], c='green'  , label='longitudinal')
    ax2.plot(recv_time,recv_lobe_emag_rad[0], c='magenta', label='radial')
    ax2.legend(framealpha=0.4)
    ax2.set_xlabel('Time (Myr)')
    ax2.set_ylabel('Lobe magnetic energy (erg)')
    fig.savefig("hardcastle14_plots/energy_test.png")