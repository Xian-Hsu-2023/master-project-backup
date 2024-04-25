import yt
import matplotlib.pyplot as plt
import numpy as np
from yt.mods import *
from yt.units import Msun, Mpc, dyne, cm, g
from yt import YTQuantity
from yt.utilities.physical_constants import kb

# file_name = "../mhd_injectiontest125_toroidal_kincoef175_ONeill_amb_inj15_ref8_0.001/"
file_name = ""
textfile_name = "../write_file_test/text_file/test129.txt"
file_number = 100
# ds = [None]*100
count = [None]*1000
t = np.zeros(100)
G = 6.67430e-8*dyne*cm**2/g**2
M = 8.5e14*Msun
r_vir = 2.440*Mpc
c = 6.81
E_k_dif = np.zeros(100)
E_th_dif = np.zeros(100)
E_th_dif_abs = np.zeros(100)
E_g_dif = np.zeros(100)
E_g_dif_abs = np.zeros(100)
E_m_dif = np.zeros(100)
E_dif =   np.zeros(100)
E_dif_abs =   np.zeros(100)

for i in range(10):
    count[i] = str('000')+str(i)
for i in range(10,100):
    count[i] = str('00')+str(i)
for i in range(100,1000):
    count[i] = str('0')+str(i)

ds_zero = yt.load(file_name+'crbub_hdf5_plt_cnt_0000')
#sp = ds_zero.sphere("c", (150, "kpc"))
sp = ds_zero.all_data()
#sp = ad.cut_region(['obj["gas", "cooling_time"].in_units("Gyr") > 3'])
E_k_initial  = sp["gas", "E_k"].sum()
E_th_initial = sp["gas", "E_th"].sum()
E_g_initial  = sp["gas", "E_g"].sum()
E_m_initial  = sp["gas", "E_m"].sum()
for i in range(0,file_number):
    ds = yt.load(file_name+'crbub_hdf5_plt_cnt_'+count[i])
    ds_base = yt.load('../mhd_injectiontest129_kincoef175_ref8_0.0/crbub_hdf5_plt_cnt_'+count[i])
    sp = ds.all_data()
    sp_base = ds_base.all_data()
    #sp = ad.cut_region(['obj["gas", "cooling_time"].in_units("Gyr") > 3'])
    t[i] = ds.current_time.in_units('Myr')
    E_k_dif[i]  = sp["gas", "E_k"].sum()-E_k_initial
    E_th_dif[i] = sp["gas", "E_th"].sum()-E_th_initial
    # E_th_dif_abs[i] = abs(E_th_dif[i])
    E_g_dif[i]  = sp["gas", "E_g"].sum()-E_g_initial
    # E_g_dif_abs[i] = abs(E_g_dif[i])
    E_m_base = sp_base["gas", "E_m"].sum()
    E_m_dif[i]  = sp["gas", "E_m"].sum()-E_m_base
    E_dif[i]    = E_k_dif[i] + E_th_dif[i] + E_g_dif[i] + E_m_dif[i]
    # E_dif_abs[i]= abs(E_dif[i])

f = open(textfile_name, "w")
f.write("k th g m all t\n")
for i in range(file_number):
    f.write("%.3e " %E_k_dif[i])
f.write("\n")
for i in range(file_number):
    f.write("%.3e " %E_th_dif[i])
f.write("\n")
for i in range(file_number):
    f.write("%.3e " %E_g_dif[i])
f.write("\n")
for i in range(file_number):
    f.write("%.3e " %E_m_dif[i])
f.write("\n")
for i in range(file_number):
    f.write("%.3e " %E_dif[i])
f.write("\n")
for i in range(file_number):
    f.write("%.3e " %t[i])
f.close()

