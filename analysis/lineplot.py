import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 12, 'font.weight' : 'bold'})
import numpy as np
import matplotlib.pyplot as plt
import yt
from yt import YTQuantity
def _faraday_rotation_y(field, data):
    coef = 0.812*3.2407792896664E-19 
    mue = 1.18 
    mp = 1.67e-24*YTQuantity(1.,"g") 
    ne = data["density"]/mue/mp
    return coef * ne * data["magnetic_field_y"] * YTQuantity(1e6, 'rad/m**2*cm**3/gauss/cm')
yt.add_field(("gas", "faraday_rotation_y"), function = _faraday_rotation_y, sampling_type="local", units="rad/m**2/cm")
KPC = 3.08567758e21 # cm

diretory = 's041_regular_r3000b8_ambLB100_ref7_0.1'
datafile = '/Data_000000'
simdir = '/lfs/data/xianhsu/gamer/bin/mhd_injection/' + diretory + datafile
UNIT_L = 1000 # kpc
ds = yt.load(simdir)
_, c = ds.find_max(("gas", "faraday_rotation_y")) # return (value, [x, y, z])
print(c*UNIT_L)
ax = 1 # y-axis
# ray = ds.ortho_ray(ax, (c[0], c[2]))
ray = ds.ortho_ray(ax, (4, 4))
srt = np.argsort(ray["index", "y"])
dy = np.array(ray["gas", "dy"][srt])
faraday_rotation_y = np.array(ray["gas", "faraday_rotation_y"][srt])
cumulative_faraday_rotation_y = np.cumsum(faraday_rotation_y*dy)
# cumulative_dy = np.cumsum(dy)
# print(np.array(ray["index", "y"][srt])*UNIT_L)
# print(cumulative_dy)
plt.subplot(311)
plt.semilogy(np.array(ray["index", "y"][srt])*UNIT_L, np.array(ray["gas", "density"][srt]), '.')
plt.ylabel(r"dens $(g\,cm^{-3})$")
plt.title('time = %.2f Myr' %ds.current_time.in_units('Myr'))
plt.subplot(312)
plt.plot(np.array(ray["index", "y"][srt])*UNIT_L, np.array(ray["gas", "magnetic_field_y"][srt]), '.')
plt.ylabel(r"B field (gauss)")
# plt.ylim(1e-7,1e-3)
plt.ylim(-1e-5, 1e-5)
plt.subplot(313)
# plt.semilogy(np.array(ray["index", "y"][srt])*UNIT_L, faraday_rotation_y, '.')
plt.semilogy(np.array(ray["index", "y"][srt])*UNIT_L, -cumulative_faraday_rotation_y, '.')
plt.ylabel(r"$FR y (rad\, m^{-2})$")
plt.ylim(1e2,1e5)
plt.xlabel("y (kpc)")
plt.savefig("LinePlot/" + diretory[0:4] + ".png")
plt.xlim(3500,4500)
plt.savefig("LinePlot/" + diretory[0:4] + "_500.png")
plt.xlim(3900,4100)
plt.savefig("LinePlot/" + diretory[0:4] + "_100.png")