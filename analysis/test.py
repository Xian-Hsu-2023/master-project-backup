import yt
from yt.utilities.physical_constants import kb, pi
# from yt.visualization.base_plot_types import get_multi_plot
from yt import YTQuantity
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import AxesGrid
import os

file_directory = ['singletest_ref1',
                  'singletest_ref2',
                  'singletest_ref3']
# plots = [None]*3
# orient = 'horizontal'
def _tCool(field, data): 
    mu = 0.61 # assuming fully ionized gas 
    mue = 1.18 
    mp = 1.67e-24*YTQuantity(1.,"g") 
    ne = data["gas", "density"]/mue/mp 
    n  = data["gas", "density"]/mu/mp 
    T  = data["gas", "temperature"].in_units('K') 
    return 1.5*n/ne**2.*kb*T**0.5/YTQuantity(3e-27, 'erg/s*cm**3/K**0.5')
yt.add_field(("gas", "cooling_time"), function = _tCool, sampling_type="local", units="Gyr")

fig = plt.figure()
plt.rcParams.update({'font.size': 10})
grid = AxesGrid(
    fig,
    (0.1, 0.1, 0.8, 0.9),
    nrows_ncols=(1, 3),
    axes_pad=0.01,
    # axes_pad=0.05,
    label_mode="L",
    share_all=True,
    cbar_location="right",
    cbar_mode="single",
    cbar_size="3%",
    cbar_pad="0%",
)
titles=['lv1','lv2','lv3']
# fig, axes, colorbars = get_multi_plot(3, 1, colorbar=orient, bw = 4)
# f_axes = [axes[0][0], axes[0][1], axes[0][2]]
for j in range(len(file_directory)):
    file_name = '../'+file_directory[j]+'/Data_000100'
    #ds = yt.load('../crbub_hdf5_plt_cnt_'+count[i])
    ds = yt.load(file_name)
    p = yt.SlicePlot(ds, 'x', 'cooling_time')
    p.set_zlim(("gas", "cooling_time"), 1, 2e2)
    p.annotate_grids()
    plot = p.plots[("gas", "cooling_time")]
    plot.figure = fig
    plot.axes = grid[j].axes
    plot.cax = grid.cbar_axes[j]
    p.render()
    grid[j].axes.set_title(titles[j])
    # slc_frb = slc.data_source.to_frb((100, "kpc"), 800)
    # slc_var = np.array(slc_frb['cooling_time'])
    # plots[j] = f_axes[j].imshow(slc_var, origin='lower', norm=LogNorm())
    # plots[j].set_cmap("gnuplot")
    # plots[j].set_clim((1,2e2))
    # plots[j].annotate_grids()
grid[0].axes.text(0.3, 0.9, 'time = %.2f Myr' %ds.current_time.in_units('Myr'), horizontalalignment='center',verticalalignment='center', transform=grid[0].axes.transAxes)
# axes[0][0].text(20,30,'time = %.2f Myr' %ds.current_time.in_units('Myr'),color='w',size=15)

try:
    fig.savefig('test.png')
except:
    print("error happens when saving file... QQ")
