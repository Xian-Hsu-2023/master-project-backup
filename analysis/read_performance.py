import matplotlib.pyplot as plt
import numpy as np

dir_name = ["efficiency_2task"
            , "efficiency_2task_1GPU"
            , "efficiency_3task"
            , "efficiency_4task"
            , "singletest_noadjbh_014_smallrinj_ref4"]

output = "efficiency.png"

fig = plt.figure(figsize=(10,7.5))
plt.rcParams.update({'font.size': 12})
ax = fig.add_subplot(111)
for i in range(len(dir_name)):
    file_name = "../" + dir_name[i] + "/Record__Performance"
    with open(file_name) as textFile:
        lines = [line.split() for line in textFile]
    data = np.array(lines[2:], dtype = float)
    time = np.array(data[:,0], dtype = float)
    perf_overall = np.array(data[:,6], dtype = float)
    # ax.plot(time, perf_overall, alpha=0.7)
    ax.scatter(time, perf_overall, alpha=0.4, s=14)
    # ax.semilogy(time, powerinj)
    
# ax.set_xlim(-0.01,1)
# ax.set_xlim(-1e-3,1e-3)
# ax.set_ylim(1e41, 1e45)
ax.set_xlabel("Time (SimUnit=250Myr, Injection Time = 0.1)")
ax.set_ylabel("Overall Performance (cells / sec)")
ax.legend(dir_name)
ax.grid()
fig.savefig(output)