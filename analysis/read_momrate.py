import matplotlib.pyplot as plt
import numpy as np

dir_name = ["singletest_noadjbh_013_smallrinj_ref3",
            "singletest_noadjbh_016_copyref3init_ref4"]

output = "momrate.png"

fig = plt.figure()
ax = fig.add_subplot(111)
for i in range(len(dir_name)):
    file_name = "../" + dir_name[i] + "/Record__Center"
    with open(file_name) as textFile:
        lines = [line.split() for line in textFile]
    data = np.array(lines[1:], dtype = float)
    time = np.array(data[:,0], dtype = float)
    injmomx = np.array(data[:,20], dtype = float)
    rinjmomx = np.zeros(len(injmomx)-1)
    for j in range(len(rinjmomx)):
        rinjmomx[j] = (injmomx[j+1] - injmomx[j]) / (time[j+1] - time[j])
    ax.loglog(time[1:], rinjmomx)

ax.set_xlabel("Time (simunit=250Myr)")
ax.set_ylabel("x-momentum injection rate (g*cm/s**2)")
ax.legend(dir_name)
ax.grid()
fig.savefig(output)