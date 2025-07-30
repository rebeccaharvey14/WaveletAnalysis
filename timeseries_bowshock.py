import sys
import datetime
import matplotlib
import numpy as np
import pandas as pd
import pyspedas
from pytplot import get_data
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits import axisartist
from mpl_toolkits.axisartist.parasite_axes import HostAxes, ParasiteAxes
from mpl_toolkits.axes_grid1 import host_subplot
from helicity_calculations import helicity
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from functions import get_variables
###########################################################################

time_head  = sys.argv[1] + ' ' + sys.argv[2]
time_tail  = sys.argv[3] + ' ' + sys.argv[4]
time_range = [time_head,time_tail]
namestr    = sys.argv[5]
probe_str  = sys.argv[6]
probe      = probe_str[4].lower()
bs_start = sys.argv[7] + ' ' + sys.argv[8]
bs_end = sys.argv[9] + ' ' + sys.argv[10]

rootDir = '/home/rharvey/Documents/Research/Wavelet-Analysis/'
dataFile = '/home/rharvey/data/' + 'data' + probe_str + namestr + '.csv'

# Read in data file
Bx, By, Bz, Vx, Vy, Vz, Np, Bmag, bmag, Vmag, Tp, Te, beta, Time, dt = get_variables(dataFile,time_range)
df = pd.read_csv(dataFile,index_col=0)
df_inre = df[['x_gse','y_gse','z_gse']]
time = np.array([t.strftime('%Y-%m-%d %H:%M:%S.%f') for t in Time])

if probe_str[:-1]=='_THM':
    pyspedas.themis.esa(probe=probe, trange=time_range, time_clip=True, varnames=['th'+probe+'_peir_en_eflux'])
    e_epoch, energies, flux = get_data('th'+probe+'_peir_en_eflux')
else: 
    pyspedas.mms.mms_load_fpi(probe=probe, trange=time_range, time_clip=True, center_measurement=True)
    e_epoch, energies, flux = get_data('mms'+probe+'_dis_energyspectr_mx_fast')
eTime = np.array([datetime.datetime.utcfromtimestamp(t) for t in e_epoch])
etime = np.array([t.strftime('%Y-%m-%d %H:%M:%S.%f') for t in eTime])

# Plot Of B-Field, Velocity, Density, Temperature, Beta, Energy
fig = plt.figure(figsize=(10,14))
fig.subplots_adjust(hspace=0)
formatter = matplotlib.dates.DateFormatter('%H:%M')
if pd.to_datetime(time_head).day == pd.to_datetime(time_tail).day:
    fig.suptitle(pd.to_datetime(str(time_head)).strftime('%d %B %Y'), fontsize=15, y=0.9)
else:
    fig.suptitle(pd.to_datetime(time_head).strftime('%d %B') + '--' + pd.to_datetime(time_tail).strftime('%d %B %Y'),y=0.9, fontsize=15)

# MAGNETIC FIELD
ax1 = host_subplot(611)
p1a = ax1.plot(Time,Bx,color='black', label='$B_X$')
p1b = ax1.plot(Time,By,color='blue', label='$B_Y$')
p1c = ax1.plot(Time,Bz,color='red', label='$B_Z$')
ax1.set_xlim([time_head,time_tail])
ax1.xaxis.set_major_formatter(formatter)
ax1.set_ylabel('(a) \n B [nT]', fontsize=15)
ax1.axes.get_xaxis().set_visible(False)
leg1 = ax1.legend(bbox_to_anchor=(1.06,0.5), loc='center right', handlelength=0, handletextpad=0, frameon=False, prop={'size': 12}, fontsize=15, labelcolor='linecolor')
ax1.axvspan(bs_start, bs_end, color='gray', alpha=0.2, lw=0)
ax1.axvline(bs_start, color='black', linewidth=0.2)
ax1.axvline(bs_end, color='black', linewidth=0.2)

# VELOCITY
ax2 = host_subplot(612)
p2a = ax2.plot(Time,Vx,color='black', label='$V_X$')
p2b = ax2.plot(Time,Vy,color='blue', label='$V_Y$')
p2c = ax2.plot(Time,Vz,color='red', label='$V_Z$')
ax2.set_xlim([time_head,time_tail])
ax2.xaxis.set_major_formatter(formatter)
ax2.set_ylabel('(b) \n V [km/s]', fontsize=15)
ax2.axes.get_xaxis().set_visible(False)
leg2 = ax2.legend(bbox_to_anchor=(1.06,0.5), loc='center right', handlelength=0, handletextpad=0, frameon=False, prop={'size': 12}, fontsize=15, labelcolor='linecolor')
ax2.axvspan(bs_start, bs_end, color='gray', alpha=0.2, lw=0)
ax2.axvline(bs_start, color='black', linewidth=0.2)
ax2.axvline(bs_end, color='black', linewidth=0.2)

# DENSITY & TEMPERATURE
ax3 = host_subplot(613)
p3 = ax3.plot(Time,Np,color='black', label='Np')
ax3a = ax3.twinx()
#p3a= ax3.plot(Time,Tp*(1e-6),color='blue', label='T')
p2a = ax3.plot(Time,Tp/11604.518,color='blue', label='T')
ax3.set_xlim([time_head,time_tail])
ax3.xaxis.set_major_formatter(formatter)
ax3.set_ylabel('(c) \n $N_p$ [cm$^{-3}$]', fontsize=15)
ax3a.set_ylabel('T [eV]', fontsize=15)
ax3.axes.get_xaxis().set_visible(False)
leg3 = ax3.legend(loc='best', handlelength=0, handletextpad=0, frameon=False, ncol=2, prop={'size': 12}, fontsize=15, labelcolor='linecolor')
ax3.axvspan(bs_start, bs_end, color='gray', alpha=0.2, lw=0)
ax3.axvline(bs_start, color='black', linewidth=0.2)
ax3.axvline(bs_end, color='black', linewidth=0.2)

# ALFVEN SPEED MAGNITUDE & FLUCTUATIONS
ax4 = host_subplot(614)
ax4.plot(Time,bmag,color='black',label='$V_A$')
ax4a = ax4.twinx()
ax4a.plot(Time,Vmag-np.nanmean(Vmag),color='blue',label=r'$\delta V$')
ax4.set_xlim([time_head,time_tail])
ax4.xaxis.set_major_formatter(formatter)
ax4.set_ylabel('(d) \n $V_A$ [km/s]', fontsize=15)
ax4a.set_ylabel(r'|V| - <V> [km/s]', fontsize=13)
ax4.axes.get_xaxis().set_visible(False)
leg4 = ax4.legend(loc='best', handlelength=0, handletextpad=0, frameon=False, ncol=2, prop={'size': 12}, fontsize=15, labelcolor='linecolor')
ax4.axvspan(bs_start, bs_end, color='gray', alpha=0.2, lw=0)
ax4.axvline(bs_start, color='black', linewidth=0.2)
ax4.axvline(bs_end, color='black', linewidth=0.2)

# PLASMA BETA
ax5 = host_subplot(615)
ax5.plot(Time,beta,color='black',label='$\\beta$')
ax5.set_xlim([time_head,time_tail])
ax5.xaxis.set_major_formatter(formatter)
ax5.set_yscale('log',base=10)
ax5.set_ylabel('(e) \n $\\beta$', fontsize=15)
ax5.axes.get_xaxis().set_visible(False)
ax5.axvspan(bs_start, bs_end, color='gray', alpha=0.2, lw=0)
ax5.axvline(bs_start, color='black', linewidth=0.2)
ax5.axvline(bs_end, color='black', linewidth=0.2)

# ION SPECTRA
ax6 = host_subplot(616)
levels = np.logspace(3,9,num=50)
cbar_ax = fig.add_axes([0.91,0.111,0.01,0.125])

#cs = ax6.contourf(eTime,flux[0], energies.T, levels=levels, locator=ticker.LogLocator(), cmap=plt.cm.nipy_spectral)
cs = ax6.contourf(eTime,flux[0][1:-8], energies.T[1:-8], levels=levels, locator=ticker.LogLocator(), cmap=plt.cm.nipy_spectral)

cbar = fig.colorbar(cs, ticks=np.logspace(3,9,num=7),format=ticker.LogFormatterSciNotation(),cax=cbar_ax)
cbar.ax.set_ylabel('Energy flux\n [eV/(cm$^2$ s sr eV)]', fontsize=15)
ax6.set_ylabel('(f) \n Ions [eV]', fontsize=15)
ax6.set_xlim([time_head,time_tail])
ax6.set_yscale('log',base=10)

ticklabel = []
for time_idx in time[np.linspace(0,time.size-1,num=6).astype(int)]:
    ticklabel.append(time_idx[11:16] + '\n %.2f\n %.2f\n %.2f' %(df_inre['x_gse'][time_idx], df_inre['y_gse'][time_idx], df_inre['z_gse'][time_idx]))

ax6.set_xticks(df.index[np.linspace(0,time.size-1,num=6).astype(int)])
ax6.set_xticklabels(ticklabel,fontsize=15)
ax6.text(-0.05,-0.54,'HH:MM\n X ($R_E$)\n Y ($R_E$)\n Z ($R_E$)',horizontalalignment='right',transform=ax6.transAxes, fontsize=15)

plt.savefig(rootDir + 'Plots/timeseries_bowshock' + namestr + '_THM' + probe.upper() + '.png', bbox_inches='tight', dpi=300)
plt.show()
plt.close()