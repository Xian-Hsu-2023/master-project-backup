import yt
from yt.utilities.physical_constants import kb
from yt import YTQuantity
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def _tCool(field, data): 
    mu = 0.61 # assuming fully ionized gas 
    mue = 1.18 
    mp = 1.67e-24*YTQuantity(1.,"g") 
    ne = data["gas", "density"]/mue/mp 
    n  = data["gas", "density"]/mu/mp 
    T  = data["gas", "temperature"].in_units('K') 
    return 1.5*n/ne**2.*kb*T**0.5/YTQuantity(3e-27, 'erg/s*cm**3/K**0.5')
yt.add_field(("gas", "cooling_time"), function = _tCool, sampling_type="local", units="Gyr") 

start_ID = int(1)
max_ID = int(100)
interval = int(1)
number = int((max_ID-start_ID)/interval+1)
dir_name = ["../s046_reg_r3kb8_LB50_2A05_ref9_0.1/",
            "../s048_reg_r3kb8_LB50_2A05_ref8_0.1/", 
            "../s050_reg_r3kb8_LB50_2A05_ref7_0.1/",
            "../s052_reg_r3kb8_LB50_2A05_ref6_0.1/"]

count = [None]*101
for i in range(101):
    count[i] = str(i).zfill(6)
# t = [None]*number
time = np.zeros(number)
length = np.zeros((len(dir_name), number))
recv_len = np.zeros_like(length)
recv_tim = np.zeros_like(time)

for j in range(len(dir_name)):
    for i in range(number):
        if (j*len(dir_name)+i) % size == rank:
            ds = yt.load(dir_name[j]+'Data_'+count[i*interval+start_ID])
            # ad = ds.all_data()
            ad = ds.sphere("c", (200,"kpc"))
            sp = ad.cut_region(['obj["gas", "cooling_time"].in_units("Gyr") > 25']) # larger than 25 Gyr is necessary
            # print(f"in time = {output1}, 20, 25, 30 Gyr cooling time obtains {output2}, {output3}, {output4}")
            # time = np.zeros(number)
            # length = np.zeros(number)
            # print(sp.argmax(("gas", "x")))
            length[j, i] = sp.argmax(("gas", "x"))[0] / YTQuantity(1, 'kpc') - 4000
            # length[j, i] = (sp.argmax(("gas", "x"))[0] -YTQuantity(500, 'kpc')) *1000
            # length[int(i/interval)] = sp.argmax(("gas", "z"))[2]/3.0857e21
            if j == 0:
                time[i] = ds.current_time.in_units('Myr')

# for j in range(len(dir_name)):
comm.Barrier()
comm.Reduce(length, recv_len)
comm.Reduce(time, recv_tim)
if rank == 0:
    # print(sp.argmax(("gas", "x")))
    print(recv_len, flush=True)
    print(recv_tim, flush=True)
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    for i in range(len(dir_name)):
        ax.plot(recv_tim,recv_len[i], alpha=0.5)
    #ax.legend(t)
    ax.set_xlim(1,250)
    ax.set_ylim(1,200)
    ax.set_xlabel('time (Myr)')
    ax.set_ylabel('jet length (kpc)')
    # ax.legend([r"$f_B=0.1\%$",r"$f_B=1\%$",r"$f_B=10\%$"])
    ax.legend(['lv9','lv8','lv7','lv6'])
    plt.title('Jet Length Evolution')
    fig.savefig('jet_length/s046.png')

# mpiexec -n 32 python estimate_jet_length_mpi.py
# mpiexec -n 6 python estimate_jet_length_mpi.py

