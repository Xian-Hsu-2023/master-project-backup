import yt
from yt.utilities.physical_constants import kb, pi
from yt.visualization.base_plot_types import get_multi_plot
from yt.mods import *
import numpy as np
import matplotlib.colorbar as cb
from matplotlib.colors import LogNorm
import sys
import os

sys.path.append(os.path.abspath("/data/xianhsu/CRBUB_GPU"))
from inputfile import *
plot_file_number = get_plot_file_numbers()
file_name_type = get_file_name_type()
#directory_name = ''
#plot_file_number = get_plot_file_numbers()
count = [None]*1000
#ds = [None]*3
# slc = [None]*4
# slc_frb = [None]*4
# slc_var = [None]*4
plots = [None]*4
orient = 'horizontal'

for i in range(1000):
    count[i] = str(i).zfill(6)
file_directory = ['singletest004_ref1',
                  'singletest003_ref2',
                  'singletest002_ref3',
                  'singletest001_ref4']

aregion = get_region()
adirection = get_direction()

def _tCool(field, data): 
    mu = 0.61 # assuming fully ionized gas 
    mue = 1.18 
    mp = 1.67e-24*YTQuantity(1.,"g") 
    ne = data["gas", "density"]/mue/mp 
    n  = data["gas", "density"]/mu/mp 
    T  = data["gas", "temperature"].in_units('K') 
    return 1.5*n/ne**2.*kb*T**0.5/YTQuantity(3e-27, 'erg/s*cm**3/K**0.5')
yt.add_field(("gas", "cooling_time"), function = _tCool, sampling_type="local", units="Gyr") 

for region in aregion:
    for direction in adirection:
        directory = 'coolingtime_'+direction+'_range'+str(region)
        command = 'mkdir '+directory
        os.system(command)
        
        for i in range(plot_file_number):    
            fig, axes, colorbars = get_multi_plot(4, 1, colorbar=orient, bw = 5)
            f_axes = [axes[0][0], axes[0][1], axes[0][2], axes[0][3]]
            for fax in f_axes:
                fax.xaxis.set_visible(False)
                fax.yaxis.set_visible(False)
            for j in range(len(file_directory)):
                file_name = '../'+file_directory[j]+'/'+file_name_type+count[i]
                #ds = yt.load('../crbub_hdf5_plt_cnt_'+count[i])
                ds = yt.load(file_name)
                slc = yt.SlicePlot(ds, direction, 'cooling_time')
                slc_frb = slc.data_source.to_frb((region, "kpc"), 800)
                slc_var = np.array(slc_frb['cooling_time'])
                plots[j] = f_axes[j].imshow(slc_var, origin='lower', norm=LogNorm())
                plots[j].set_cmap("gnuplot")
                plots[j].set_clim((1,2e2))
            # plots[0].set_clim((1e-27,2e-25))
            # plots[1].set_clim((3e-11,5e-9))
            # plots[2].set_clim((3e7,1e8))
            # plots[0].set_cmap("gnuplot")
            # plots[1].set_cmap("gnuplot")
            # plots[2].set_cmap("coolwarm")
            axes[0][0].text(20,30,'time = %.2f Myr' %ds.current_time.in_units('Myr'),color='w',size=15)
                #axes[0][0].text(5,760,'density',color='w',size=15)
                #axes[0][1].text(5,760,'pressure',color='w',size=15)
                #axes[0][2].text(5,760,'temperature',color='w',size=15)
            titles=['lv1','lv2','lv3','lv4']
            for p, cax, t in zip(plots[0:5], colorbars, titles):
                cbar = fig.colorbar(p, cax=cax, orientation=orient)
                cbar.set_label(t)
            try:
                fig.savefig(directory+'/'+count[i]+'.png')
            except:
                print("error happens when saving file... QQ")
            del fig, axes, colorbars


#axes[0][0].text(5,760,'20% magnetic field energy',color='w',size=15)
#axes[0][1].text(5,760,'50% magnetic field energy',color='w',size=15)
#axes[0][2].text(5,760,'80% magnetic field energy',color='w',size=15)
#fig.suptitle('Density with different magnetic field energy ratio',color='w',fontsize=20)
#titles=[r'$\mathrm{B0%}\ (\mathrm{g\ cm^{-3}})$', r'$\mathrm{B50%}\ (\mathrm{g\ cm^{-3}})$', r'$\mathrm{B90%}\ (\mathrm{g\ cm^{-3}})$']
#titles=['Density at 20%B','Density at 50%B','Density at 80%B']
#for p, cax, t in zip(plots[0:12:4], colorbars, titles):
#    cbar = fig.colorbar(p, cax=cax, orientation=orient)
#    cbar.set_label(t)
#fig.savefig('001.png')
