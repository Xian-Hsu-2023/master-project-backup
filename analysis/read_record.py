import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import re
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 12, 'font.weight' : 'bold'})

directory = "../s057_t250_r3kb8_LB50_ref7_0.0"
textfile_name = directory + "/Record__Conservation"
base_directory = "../s051_r3kb8_LB50_ref7"
base_textfile_name = base_directory + "/Record__Conservation"
# print(directory[2:8])
# output = 'evolution_plot/' + 's021_record.png' # TODO: to be modified
output = 'evolution_plot' + directory[2:8] + 'record.png' 
# output = 'evolution_plot/' + 'past_record_gas.png'
# if 'base_directory' in globals(): # test: success
#     print('1')

# f_b = 0.01
f_b = float(re.findall('[\d.]+',directory)[-1])
jet_power = 5e45
jet_duration = 250
MYR = 86400 * 365.25 * 1e6
UNIT_L = 3.08567758149e24
UNIT_M = 1.9885e47
UNIT_T = 7.8892313e15
UNIT_E = UNIT_M * UNIT_L**2 / UNIT_T**2
# print(UNIT_E)
bin_number = int(101)
# plot variables
xmin, xmax = 0, 150
ymin, ymax = 1e56, 1e63

# note! (simulation Record__Conservation with mhd)
# total energy     => 27 + 51 = 66 (hydro+par=all)
# kinetic energy   => 15 + 45 (hydro+par)
# thermal energy   => 18
# potential energy => 21 + 48 (hydro+par)
# magnetic energy  => 24

with open(textfile_name) as textFile:
    lines = [line.split() for line in textFile]
data          = np.array(lines[27:], dtype = float)

time = np.array(data[:, 0]) * UNIT_T / MYR
ekin_gas = (np.array(data[:,14]) - np.array(data[0,14])) * UNIT_E
ekin_par = (np.array(data[:,44]) - np.array(data[0,44])) * UNIT_E
ethe_gas = (np.array(data[:,17]) - np.array(data[0,17])) * UNIT_E
epot_gas = (np.array(data[:,20]) - np.array(data[0,20])) * UNIT_E
epot_par = (np.array(data[:,47]) - np.array(data[0,47])) * UNIT_E
emag_gas = (np.array(data[:,23]) - np.array(data[0,23])) * UNIT_E
etot_gas = (np.array(data[:,26]) - np.array(data[0,26])) * UNIT_E
etot_par = (np.array(data[:,50]) - np.array(data[0,50])) * UNIT_E
etot_all = (np.array(data[:,65]) - np.array(data[0,65])) * UNIT_E
legend = ['kinetic gas', 'kinetic par', 'thermal gas', 
          'gravitational gas', 'gravitational par', 'magnetic gas', 
          'total gas', 'total par', 'total all', 'magnetic injection', 'total injection']

# TODO: determine the interpolated value (in order to do subtraction between two simulations)

file_time = np.linspace(0, 1, num=bin_number) * UNIT_T / MYR
file_ekin_gas = np.zeros(bin_number)
file_ekin_par = np.zeros(bin_number)
file_ethe_gas = np.zeros(bin_number)
file_epot_gas = np.zeros(bin_number)
file_epot_par = np.zeros(bin_number)
file_emag_gas = np.zeros(bin_number)
file_etot_gas = np.zeros(bin_number)
file_etot_par = np.zeros(bin_number)
file_etot_all = np.zeros(bin_number)
count = int(0)
for i in range(len(time)-1):
    fin = True
    while fin:
        if (file_time[count]>=time[i] and file_time[count]<=time[i+1]):
            r2 = (file_time[count]-time[i]) / (time[i]-time[i+1])
            r1 = 1 - r2
            file_ekin_gas[count] = ekin_gas[i]*r1 + ekin_gas[i+1]*r2
            file_ekin_par[count] = ekin_par[i]*r1 + ekin_par[i+1]*r2
            file_ethe_gas[count] = ethe_gas[i]*r1 + ethe_gas[i+1]*r2
            file_epot_gas[count] = epot_gas[i]*r1 + epot_gas[i+1]*r2
            file_epot_par[count] = epot_par[i]*r1 + epot_par[i+1]*r2
            file_emag_gas[count] = emag_gas[i]*r1 + emag_gas[i+1]*r2
            file_etot_gas[count] = etot_gas[i]*r1 + etot_gas[i+1]*r2
            file_etot_par[count] = etot_par[i]*r1 + etot_par[i+1]*r2
            file_etot_all[count] = etot_all[i]*r1 + etot_all[i+1]*r2
            count += 1
            if count==bin_number:
                break
        else:
            fin = False

# for information in [file_time, file_ekin_gas, file_ekin_par, file_ethe_gas, file_epot_gas, 
#                     file_epot_par, file_emag_gas, file_etot_gas, file_etot_par, file_etot_all]:
#     print(information)

with open(base_textfile_name) as textFile:
    lines = [line.split() for line in textFile]
data          = np.array(lines[27:], dtype = float)

time = np.array(data[:, 0]) * UNIT_T / MYR
ekin_gas = (np.array(data[:,14]) - np.array(data[0,14])) * UNIT_E
ekin_par = (np.array(data[:,44]) - np.array(data[0,44])) * UNIT_E
ethe_gas = (np.array(data[:,17]) - np.array(data[0,17])) * UNIT_E
epot_gas = (np.array(data[:,20]) - np.array(data[0,20])) * UNIT_E
epot_par = (np.array(data[:,47]) - np.array(data[0,47])) * UNIT_E
emag_gas = (np.array(data[:,23]) - np.array(data[0,23])) * UNIT_E
etot_gas = (np.array(data[:,26]) - np.array(data[0,26])) * UNIT_E
etot_par = (np.array(data[:,50]) - np.array(data[0,50])) * UNIT_E
etot_all = (np.array(data[:,65]) - np.array(data[0,65])) * UNIT_E

# file_time = np.linspace(0, 1, num=bin_number) * UNIT_T / MYR
base_ekin_gas = np.zeros(bin_number)
base_ekin_par = np.zeros(bin_number)
base_ethe_gas = np.zeros(bin_number)
base_epot_gas = np.zeros(bin_number)
base_epot_par = np.zeros(bin_number)
base_emag_gas = np.zeros(bin_number)
base_etot_gas = np.zeros(bin_number)
base_etot_par = np.zeros(bin_number)
base_etot_all = np.zeros(bin_number)
count = int(0)
for i in range(len(time)-1):
    fin = True
    while fin:
        if (file_time[count]>=time[i] and file_time[count]<=time[i+1]):
            r2 = (file_time[count]-time[i]) / (time[i]-time[i+1])
            r1 = 1 - r2
            base_ekin_gas[count] = ekin_gas[i]*r1 + ekin_gas[i+1]*r2
            base_ekin_par[count] = ekin_par[i]*r1 + ekin_par[i+1]*r2
            base_ethe_gas[count] = ethe_gas[i]*r1 + ethe_gas[i+1]*r2
            base_epot_gas[count] = epot_gas[i]*r1 + epot_gas[i+1]*r2
            base_epot_par[count] = epot_par[i]*r1 + epot_par[i+1]*r2
            base_emag_gas[count] = emag_gas[i]*r1 + emag_gas[i+1]*r2
            base_etot_gas[count] = etot_gas[i]*r1 + etot_gas[i+1]*r2
            base_etot_par[count] = etot_par[i]*r1 + etot_par[i+1]*r2
            base_etot_all[count] = etot_all[i]*r1 + etot_all[i+1]*r2
            count += 1
            if count==bin_number:
                break
        else:
            fin = False

diff_ekin_gas = file_ekin_gas
diff_ekin_par = file_ekin_par
diff_ethe_gas = file_ethe_gas
diff_epot_gas = file_epot_gas
diff_epot_par = file_epot_par
diff_emag_gas = file_emag_gas
diff_etot_gas = file_etot_gas
diff_etot_par = file_etot_par
diff_etot_all = file_etot_all

diff_ekin_gas = file_ekin_gas - base_ekin_gas
diff_ekin_par = file_ekin_par - base_ekin_par
diff_ethe_gas = file_ethe_gas - base_ethe_gas
diff_epot_gas = file_epot_gas - base_epot_gas
diff_epot_par = file_epot_par - base_epot_par
diff_emag_gas = file_emag_gas - base_emag_gas
diff_etot_gas = file_etot_gas - base_etot_gas
diff_etot_par = file_etot_par - base_etot_par
diff_etot_all = file_etot_all - base_etot_all
diff_plot = [diff_ekin_gas, diff_ekin_par, diff_ethe_gas, 
             diff_epot_gas, diff_epot_par, diff_emag_gas, 
             diff_etot_gas, diff_etot_par, diff_etot_all]

# print out diff_plot elements tests
# print(len(diff_plot))
# print(diff_plot)
# print(diff_plot[2])
# print(diff_plot[2][1])
# print(diff_plot[2][2])
# print(diff_plot[2][3])


fig = plt.figure()
ax = fig.add_subplot(111)
for e in diff_plot:
# for e in [ekin, ethe, epot, etot]:
    ax.semilogy(file_time, e)
# print(len(file_time), len(data[:,0]))
inj_tot = np.array([jet_power*min(jet_duration,file_time[i])*MYR for i in range(len(data[:,0]))])
inj_mag = f_b * inj_tot
# ax.plot(data[:,0], inj_mag, '--')
# ax.plot(data[:,0], inj_tot, '--')
ax.semilogy(file_time, inj_mag, '--')
ax.semilogy(file_time, inj_tot, '--')
lg = ax.legend(legend, loc='lower right')
lg.get_frame().set_alpha(0.5)
ax.set_xlabel('time (Myr)')
ax.set_ylabel(r'$E(t)-E(t_0)-E_{base}(t)$ (erg)')
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
# ax.set_ylim(1e56, 1e63)
ax.grid()
ax.set_title('simulation box energy evolution with particle')
fig.savefig(output)

del fig, ax

fig = plt.figure()
ax = fig.add_subplot(111)
negative_ekin_gas = np.zeros(bin_number)
negative_ekin_par = np.zeros(bin_number)
negative_ethe_gas = np.zeros(bin_number)
negative_epot_gas = np.zeros(bin_number)
negative_epot_par = np.zeros(bin_number)
negative_emag_gas = np.zeros(bin_number)
negative_etot_gas = np.zeros(bin_number)
negative_etot_par = np.zeros(bin_number)
negative_etot_all = np.zeros(bin_number)
negative_plot = [negative_ekin_gas, negative_ekin_par, negative_ethe_gas, 
                 negative_epot_gas, negative_epot_par, negative_emag_gas, 
                 negative_etot_gas, negative_etot_par, negative_etot_all]
for i in range(bin_number):
    for j in range(len(diff_plot)):
        if diff_plot[j][i] < 0 :
            negative_plot[j][i] = -diff_plot[j][i]

for e in negative_plot:
# for e in [ekin, ethe, epot, etot]:
    ax.semilogy(file_time, e)
inj_tot = np.array([jet_power*min(jet_duration,file_time[i])*MYR for i in range(len(data[:,0]))])
inj_mag = f_b * inj_tot
# ax.plot(data[:,0], inj_mag, '--')
# ax.plot(data[:,0], inj_tot, '--')
ax.semilogy(file_time, inj_mag, '--')
ax.semilogy(file_time, inj_tot, '--')
lg = ax.legend(legend, loc='lower right')
lg.get_frame().set_alpha(0.5)
ax.set_xlabel('time (Myr)')
ax.set_ylabel(r'$E(t)-E(t_0)-E_{base}(t)$ (erg)')
# ax.set_ylabel(r'$E(t)-E_0(t)$ (erg)')
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
# ax.set_ylim(1e56, 1e62)
# ax.set_ylim(1e56, 1e63)
ax.grid()
ax.set_title('simulation box energy evolution with particle (negative part)')
# output = 'evolution_plot/' + 's021_record_negative.png'
output = 'evolution_plot' + directory[2:8] + 'record_negative.png' 
fig.savefig(output)