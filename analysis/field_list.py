import yt
# from yt.mods import *

#ds = yt.load('../mhd_usecool_ref50_B0.5/crbub_hdf5_plt_cnt_0010')
# ds = yt.load('/lfs/data/xianhsu/gamer/bin/weinberger/old_simulations/efficiency_2task/Data_000000')
ds = yt.load('/lfs/data/xianhsu/gamer/bin/mhd_injection/s001_firsttest_ref5_0.001/Data_000000')
print(ds.derived_field_list) 
