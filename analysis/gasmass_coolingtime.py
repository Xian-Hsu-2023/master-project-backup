import yt
from yt.utilities.physical_constants import kb
from yt import YTQuantity
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

def _tCool(field, data): 
    mu = 0.61 # assuming fully ionized gas 
    mue = 1.18 
    mp = 1.67e-24*YTQuantity(1.,"g") 
    ne = data["gas", "density"]/mue/mp 
    n  = data["gas", "density"]/mu/mp 
    T  = data["gas", "temperature"].in_units('K') 
    return 1.5*n/ne**2.*kb*T**0.5/YTQuantity(3e-27, 'erg/s*cm**3/K**0.5')
yt.add_field(("gas", "cooling_time"), function = _tCool, sampling_type="local", units="Gyr")

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['text.usetex'] = 'true'
plt.rcParams.update({'font.size': 16})
fig = plt.figure()
ax = fig.add_subplot(111)

dir = np.array(['../s020_compare_ref3',
                '../s021_compare_ref4',
                '../s022_compare_ref5',
                '../s023_compare_ref6'])
label = np.array(['lv3', 'lv4', 'lv5', 'lv6'])
dataname = "/Data_000020"

for i in range(len(dir)):
    ds = yt.load(dir[i]+dataname)
    sp = ds.sphere("c", (200,"kpc"))
    prof = yt.create_profile(
        sp,
        "cooling_time",
        ("gas", "mass"),
        # units={"radius": "kpc"},
        extrema={"cooling_time": ((1, "Myr"), (10000.0, "Myr"))},
        # weight_field=("gas", "mass")
        )
    colt = prof.x
    mass = prof["gas", "mass"]
    ax.semilogy(colt, mass)

ax.legend(label)
ax.set_xlabel('cooling time (Gyr)')
ax.set_ylabel('mass (g)')
ax.set_title('Gas Cooling Time Histogram')
# ax.set_xlim()
fig.savefig("ColdGasMass.png")