import yt
from yt.utilities.physical_constants import kb, pi
from yt.mods import *
import matplotlib.pyplot as plt
import numpy as np

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
max_ID = int(8)
interval = int(1)
number = int((max_ID)/interval+1)

count = [None]*5000
t = [None]*number
time = np.zeros(number)
length = np.zeros(number)
count[0] = '000000'
for i in range(1,10):
    count[i] = str('00000')+str(i)
for i in range(10,100):
    count[i] = str('0000')+str(i)
for i in range(100,1000):
    count[i] = str('000')+str(i)

fig = plt.figure()
ax  = fig.add_subplot(111)
record_ID = int(2)
for dir_name in ["../s020_compare_ref3/",
                 "../s021_compare_ref4/", 
                 "../s022_compare_ref5/",
                 "../s023_compare_ref6/",
                 "../s024_compare_ref7/"]:
    for i in range(start_ID,max_ID+1,interval):
        ds = yt.load(dir_name+'Data_'+count[i])
        # ad = ds.all_data()
        ad = ds.sphere("c", (200,"kpc"))
        sp = ad.cut_region(['obj["gas", "cooling_time"].in_units("Gyr") > 10']) # larger than 25 Gyr is necessary
        # sp1 = ad.cut_region(['obj["gas", "cooling_time"].in_units("Gyr") > 10']) # larger than 25 Gyr is necessary
        # sp2 = ad.cut_region(['obj["gas", "cooling_time"].in_units("Gyr") > 20']) # larger than 25 Gyr is necessary
        # sp3 = ad.cut_region(['obj["gas", "cooling_time"].in_units("Gyr") > 30']) # larger than 25 Gyr is necessary
        # output1 = ds.current_time.in_units('Myr')
        # output2 = sp1.argmax(("gas", "x"))[2]
        # output3 = sp2.argmax(("gas", "x"))[2]
        # output4 = sp3.argmax(("gas", "x"))[2]
        # print(f"in time = {output1}, 20, 25, 30 Gyr cooling time obtains {output2}, {output3}, {output4}")
        time = np.zeros(number)
        length = np.zeros(number)
        length[int(i/interval)] = sp.argmax(("gas", "x"))[2]
        # length[int(i/interval)] = sp.argmax(("gas", "z"))[2]/3.0857e21
        time[int(i/interval)] = ds.current_time.in_units('Myr')
    ax.loglog(time,length, alpha=0.5)
#ax.legend(t)
ax.set_xlim(1,250)
# ax.set_ylim(1e1,2e2)
# ax.set_xlabel('time (Myr)')
# ax.set_ylabel('jet length (kpc)')
# ax.legend([r"$f_B=0.1\%$",r"$f_B=1\%$",r"$f_B=10\%$"])
ax.legend(['lv3','lv2','lv3','lv4'])
plt.title('Lobe Length Evolution')
fig.savefig('jet_length.png')

