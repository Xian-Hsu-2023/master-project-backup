import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 12, 'font.weight' : 'bold'})

textfile_name = 'text_file/' + 'Record__s005_CompareToFLASH_ref4_0.001.txt'
output = 'evolution_plot/' + 's005_box_evolution.png'
legend = ['kinetic', 'thermal', 'gravitational (negative)', 'magnetic', 
          'total', 'magnetic injection', 'total injection']
f_b = 0.1
jet_duration = 10
MYR = 86400 * 365.25 * 1e6

with open(textfile_name) as textFile:
    lines = [line.split() for line in textFile]
data          = np.array(lines[2:], dtype = float)
# time kinetic thermal gravitational magnetic total

fig = plt.figure()
ax = fig.add_subplot(111)
for i in range(5): # plot basic data
    ax.semilogy(data[:,0], abs(data[:,i+1]))
inj_tot = np.array([5e45*min(jet_duration,data[i,0])*MYR for i in range(len(data[:,0]))])
inj_mag = f_b * inj_tot
ax.semilogy(data[:,0], inj_mag, '--')
ax.semilogy(data[:,0], inj_tot, '--')
lg = ax.legend(legend, loc='lower right')
lg.get_frame().set_alpha(0.5)
ax.set_xlabel('time (Myr)')
ax.set_ylabel('energy (erg)')
ax.set_xlim(0, 250)
ax.set_ylim(1e56, 1e62)
ax.grid()
ax.set_title('simulation box energy evolution')
fig.savefig(output)
