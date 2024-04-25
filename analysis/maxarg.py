import numpy as np
import yt

diretory = 's040_regular_r3000b8_amb_ref8_0.1'
datafile = '/Data_000005'
simdir = '/lfs/data/xianhsu/gamer/bin/mhd_injection/' + diretory + datafile

ds = yt.load(simdir)
sp = ds.sphere('c', (1000, "kpc"))
value, location = ds.find_max(('gas', 'sound_speed'))
print(value, location)
value, location = ds.find_max(('gas', 'temperature'))
print(value, location)
value, location = ds.find_max(('gas', 'alfven_speed'))
print(value, location)
value, location = ds.find_max(('gas', 'velocity_magnitude'))
print(value, location)
