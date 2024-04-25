import matplotlib.pyplot as plt
import numpy as np

dir_name = ["reg_ef1_ref4"
            , "reg_ef10_ref4"
            , "reg_ef100_ref4"]

output = "jet_power.png"

fig = plt.figure()
ax = fig.add_subplot(111)
for i in range(len(dir_name)):
    file_name = "../" + dir_name[i] + "/Record__Center"
    with open(file_name) as textFile:
        lines = [line.split() for line in textFile]
    data = np.array(lines[2:], dtype = float)
    time = np.array(data[:,0], dtype = float)
    powerinj = np.array(data[:,29], dtype = float)
    ax.loglog(time, powerinj)
    # ax.semilogy(time, powerinj)
    
# ax.set_xlim(-0.01,1)
# ax.set_xlim(-1e-3,1e-3)
# ax.set_ylim(1e41, 1e45)
ax.set_xlabel("Time (simunit=2Gyr)")
ax.set_ylabel("Power (erg/s)")
ax.legend(dir_name)
ax.grid()
fig.savefig(output)