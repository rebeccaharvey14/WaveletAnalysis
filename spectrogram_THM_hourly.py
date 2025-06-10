import sys
import datetime
import matplotlib
import numpy as np
import pandas as pd
import pyspedas
from pyspedas.themis.config import CONFIG
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
###########################################################################
def nan_checker(array):
	idx_nans = np.where(np.isnan(array)==True)[0]
	return idx_nans.size/array.size

def gap_checker(array,idx1,idx2,dt):
	if idx2 > Time.size:
		idx2 = Time.size-1
	elif (Time[idx2]-Time[idx1]) > datetime.timedelta(seconds=3600*3):
		while (Time[idx2]-Time[idx1]) > datetime.timedelta(seconds=800*dt*3):
			idx2 -= 1
	return idx1,idx2
###########################################################################

rootDir = '/home/rharvey/Documents/Research/Wavelet-Analysis/'
dataDir = '/home/rharvey/data/preprocessed/'
datafile = dataDir + 'DataFrame_THM' + probe.upper() + namestr + '.csv'
eventDir = rootDir + 'themis_events/'

# READ IN DATA
time_head = sys.argv[1] + ' ' + sys.argv[2]
time_tail = sys.argv[3] + ' ' + sys.argv[4]
time_range = [time_head,time_tail]
namestr = sys.argv[5]
probe = sys.argv[6]
year = sys.argv[7]

# GET EPHEMERIS DATA
coord = 'gse'
plane = 'xy'
km_in_re = 6371.2

pyspedas.themis.state(probe=probe, trange=time_range, time_clip=True, varnames=['th'+probe+'_pos_gse'])
[eph_epoch, pos, components] = get_data('th'+probe+'_pos_gse') #km
x = pos[:,0]/km_in_re
y = pos[:,1]/km_in_re
z = pos[:,2]/km_in_re

orbit_Time = np.array([datetime.datetime.utcfromtimestamp(t) for t in eph_epoch]) # dt = 1 min
orbitTime = np.array([t.strftime('%Y-%m-%d %H:%M:%S.%f') for t in orbit_Time])
titlestr = np.min(orbit_Time).strftime('%B %d %Y') + ' -- ' + np.max(orbit_Time).strftime('%B %d %Y')
df_ephem = pd.DataFrame(np.array([x,y,z]).T, index=orbitTime, columns = ['x_gse','y_gse','z_gse'])

# READ IN DATA
df = pd.read_csv(datafile,index_col=0)
selected_index = [(df.index >= time_head) & (df.index <= time_tail)][0]
df = df[selected_index]

Np = df.Np.values
Bmag = df.Bmag.values
Tp = df.Tp.values
Te = df.Te.values
Time = pd.to_datetime(df.index)
Vmag = df.Vmag.values
beta = df.beta.values
Br = -df.Bx.values
Bt = -df.By.values
Bn = df.Bz.values
Vr = -df.Vx.values
Vt = -df.Vy.values
Vn = df.Vz.values
bx = df.bx.values
by = df.by.values
bz = df.bz.values
bmag = df.bmag.values
dt = Time[1]-Time[0]
dt = dt.seconds + dt.microseconds*1e-6

Bx = df.Bx.values
By = df.By.values
Bz = df.Bz.values
Vx = df.Vx.values
Vy = df.Vy.values
Vz = df.Vz.values 

epoch = ((Time-pd.Timestamp('1970-01-01')) // pd.Timedelta('1s')).values
Time = np.array([datetime.datetime.utcfromtimestamp(t) for t in epoch])
time = np.array([t.strftime('%Y-%m-%d %H:%M:%S.%f') for t in Time])
df_inre = df_ephem.reindex(index=time,method='ffill')

pyspedas.themis.esa(probe=probe, trange=time_range, time_clip=True, varnames=['th'+probe+'_peir_en_eflux'])
e_epoch, energies, flux = get_data('th'+probe+'_peir_en_eflux')
eTime = np.array([datetime.datetime.utcfromtimestamp(t) for t in e_epoch])
etime = np.array([t.strftime('%Y-%m-%d %H:%M:%S.%f') for t in eTime])

origcd = os.getcwd()
# PLOT OF B-FIELD, VELOCITY, DENSITY, TEMPERATURE, BETA, ENERGY
fig = plt.figure(figsize=(10,14))
fig.subplots_adjust(hspace=0)
formatter = matplotlib.dates.DateFormatter('%H:%M')
if pd.to_datetime(time_head).day == pd.to_datetime(time_tail).day:
    fig.suptitle(pd.to_datetime(str(time_head)).strftime('%d %B %Y'), fontsize=15, y=0.9)
else:
    fig.suptitle(pd.to_datetime(time_head).strftime('%d %B') + '--' +pd.to_datetime(time_tail).strftime('%d %B %Y'),y=0.9, fontsize=15)
    
# MAGNETIC FIELD
ax4 = host_subplot(611)
p41 = ax4.plot(Time,Bx,color='black', label='$B_X$')
p42 = ax4.plot(Time,By,color='blue', label='$B_Y$')
p43 = ax4.plot(Time,Bz,color='red', label='$B_Z$')
ax4.set_xlim([np.nanmin(Time),np.nanmax(Time)])
ax4.xaxis.set_major_formatter(formatter)
ax4.set_ylabel('(a) \n B [nT]', fontsize=15)
ax4.axes.get_xaxis().set_visible(False)
leg4 = ax4.legend(bbox_to_anchor=(1.06,0.5), loc='center right', handlelength=0, handletextpad=0, frameon=False, prop={'size': 12}, fontsize=15)
colors=[p41[0].get_color(),p42[0].get_color(),p43[0].get_color()]
for color,text in zip(colors,leg4.get_texts()):
	text.set_color(color)

# VELOCITY
ax5 = host_subplot(612)
p51 = ax5.plot(Time,Vx,color='black', label='$V_X$')
p52 = ax5.plot(Time,Vy,color='blue', label='$V_Y$')
p53 = ax5.plot(Time,Vz,color='red', label='$V_Z$')
ax5.set_xlim([np.nanmin(Time),np.nanmax(Time)])
ax5.xaxis.set_major_formatter(formatter)
ax5.set_ylabel('(b) \n V [km/s]', fontsize=15)
ax5.axes.get_xaxis().set_visible(False)
leg5 = ax5.legend(bbox_to_anchor=(1.06,0.5), loc='center right', handlelength=0, handletextpad=0, frameon=False, prop={'size': 12}, fontsize=15)
colors=[p51[0].get_color(),p52[0].get_color(),p53[0].get_color()]
for color,text in zip(colors,leg5.get_texts()):
	text.set_color(color)

# DENSITY & TEMPERATURE
ax6 = host_subplot(613)
p61 = ax6.plot(Time,Np,color='black', label='Np')
ax7 = ax6.twinx()
#p62 = ax7.plot(Time,Tp*(1e-6),color='blue', label='T')
p62 = ax7.plot(Time,Tp/11604.518,color='blue', label='T')
ax6.set_xlim([np.nanmin(Time),np.nanmax(Time)])
ax6.xaxis.set_major_formatter(formatter)
ax6.set_ylabel('(c) \n $N_p$ [cm^-3]', fontsize=15)
ax7.set_ylabel('T [eV]', fontsize=15)
#ax7.set_yscale('log',base=10)
ax6.axes.get_xaxis().set_visible(False)
leg6 = ax6.legend(loc='best', handlelength=0, handletextpad=0, frameon=False, ncol=2, prop={'size': 12}, fontsize=15) #bbox_to_anchor=(0.9,0.64),
colors=[p61[0].get_color(),p62[0].get_color()]
for color,text in zip(colors,leg6.get_texts()):
	text.set_color(color)

# ALFVEN SPEED MAGNITUDE & FLUCTUATIONS
ax9 = host_subplot(614)
p91 = ax9.plot(Time,bmag,color='black',label='$V_A$')
ax10 = ax9.twinx()
p92 = ax10.plot(Time,Vmag-np.nanmean(Vmag),color='blue',label='$\delta V$')
ax9.set_xlim([np.nanmin(Time),np.nanmax(Time)])
ax9.xaxis.set_major_formatter(formatter)
ax9.set_ylabel('(d) \n $V_A$ [km/s]', fontsize=15)
Vstring = '|V| - {0:.0f} [km/s]'.format(np.nanmean(Vmag))
ax10.set_ylabel(r'|V| - <V> [km/s]', fontsize=13) #'V [km/s]'
ax9.axes.get_xaxis().set_visible(False)
leg9 = ax9.legend(loc='best', handlelength=0, handletextpad=0, frameon=False, ncol=2, prop={'size': 12}, fontsize=15, labelcolor='linecolor')
# colors=[p91[0].get_color(),p92[0].get_color()]
# for color,text in zip(colors,leg9.get_texts()):
# 	text.set_color(color)

# PLASMA BETA
ax11 = host_subplot(615)
p111 = ax11.plot(Time,beta,color='black',label='$\\beta$')
ax11.set_xlim([np.nanmin(Time),np.nanmax(Time)])
ax11.xaxis.set_major_formatter(formatter)
ax11.set_yscale('log',base=10)
ax11.set_ylabel('(e) \n $\\beta$', fontsize=15)
ax11.axes.get_xaxis().set_visible(False)

# ION SPECTRA
ax12 = host_subplot(616)
levels = np.logspace(3,9,num=50)
cbar_ax1 = fig.add_axes([0.91,0.111,0.01,0.125]) #[0.91,0.111,0.01,0.11]
#if int(year) < 2010:
#	cs1 = ax12.contourf(eTime,flux[0][1:flux.shape[0]], energies.T[1:flux.shape[0]],levels=levels, locator=ticker.LogLocator(), cmap=plt.cm.nipy_spectral)
#elif year=='2022' and probe=='b':
#        cs1 = ax12.contourf(eTime,flux[0][1:-16], energies.T[1:-16],levels=levels, locator=ticker.LogLocator(), cmap=plt.cm.nipy_spectral)
#else:
cs1 = ax12.contourf(eTime,flux[0][1:-8], energies.T[1:-8],levels=levels, locator=ticker.LogLocator(), cmap=plt.cm.nipy_spectral)
cbar1 = fig.colorbar(cs1, ticks=np.logspace(3,9,num=7),format=ticker.LogFormatterSciNotation(),cax=cbar_ax1)
cbar1.ax.set_ylabel('Energy flux\n [eV/(cm$^2$ s sr eV)]', fontsize=15)
ax12.set_ylabel('(f) \n Ions [eV]', fontsize=15)
ax12.set_xlim([np.nanmin(Time),np.nanmax(Time)])
ax12.set_yscale('log',base=10)

ticklabel = []
for time_idx in time[np.linspace(0,time.size-1,num=6).astype(int)]:
	ticklabel.append(time_idx[11:16] + '\n %.2f\n %.2f\n %.2f' %(df_inre['x_gse'][time_idx], df_inre['y_gse'][time_idx], df_inre['z_gse'][time_idx]))

ax12.set_xticks(df.index[np.linspace(0,time.size-1,num=6).astype(int)])
ax12.set_xticklabels(ticklabel,fontsize=15)
ax12.text(-0.05,-0.54,'HH:MM\n X ($R_E$)\n Y ($R_E$)\n Z ($R_E$)',horizontalalignment='right',transform=ax12.transAxes, fontsize=15)

plt.savefig(rootDir + '/Plots/timeseries_' + pd.to_datetime(str(time_head)).strftime('%d%m%Y') + '_THM' + probe.upper() + '.png', bbox_inches='tight', dpi=300)
#plt.show()
plt.close()
