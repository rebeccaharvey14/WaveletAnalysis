import os
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

# READ IN DATA
namestr = sys.argv[1]
probe = sys.argv[2]
rootDir = '/home/rharvey/Documents/Research/Wavelet-Analysis/'
data_dir='/home/rharvey/data/preprocessed/'

datafile = dataDir + 'DataFrame_MMS1' + namestr + '.csv'
eventFile = rootDir + 'mms_events/' + namestr[:-1] + '_MMS' + probe + '.csv'

# GET EPHEMERIS DATA
coord = 'gse'
plane = 'xy'
km_in_re = 6371.2

pyspedas.mms.mms_load_mec(probe=probe, trange=time_range, time_clip=True, varformat='*_r_' + coord)
[eph_epoch, pos] = get_data('mms1_mec_r_' + coord) 
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
dj = 0.125
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

epoch = ((Time-pd.Timestamp('1970-01-01')) // pd.Timedelta('1s')).values
Time = np.array([datetime.datetime.utcfromtimestamp(t) for t in epoch])
time = np.array([t.strftime('%Y-%m-%d %H:%M:%S.%f') for t in Time])
df_inre = df_ephem.reindex(index=time,method='ffill')
#bmag = np.array([np.sqrt(bx[i]**2 + by[i]**2 + bz[i]**2) for i in range(epoch.size)])
Bx = df.Bx.values
By = df.By.values
Bz = df.Bz.values
Vx = df.Vx.values
Vy = df.Vy.values
Vz = df.Vz.values

idx1 = 0
idx2 = idx1 + 800*3
while idx1 < Time.size-3*800 and idx2 < Time.size:

	idx1,idx2 = gap_checker(Time,idx1,idx2,dt)
	time_head = Time[idx1]
	time_tail = Time[idx2]
	print(Time[idx1],Time[idx2])

	# CHECK FOR GAPS
	if nan_checker(Br[idx1:idx2]) < 0.05 and nan_checker(Vr[idx1:idx2]) < 0.05 and nan_checker(Np[idx1:idx2]) < 0.05:

		# CALCULATE MAGNETIC HELICITY, CROSS HELICITY, AND RESIDUAL ENERGY
		scale, coi, sig_m, sig_c, sig_r = helicity(Br[idx1:idx2], Bt[idx1:idx2], Bn[idx1:idx2], Vr[idx1:idx2], Vt[idx1:idx2], Vn[idx1:idx2], Np[idx1:idx2], dt, dj)

		# PLOT OF B-FIELD AND NORMALIZED MAGNETIC HELICITY, CROSS HELICITY, & RESIDUAL ENERGY
		fig = plt.figure(figsize=(10,14))
		fig.subplots_adjust(hspace=0)
		formatter = matplotlib.dates.DateFormatter('%H:%M')
		if pd.to_datetime(time_head).day == pd.to_datetime(time_tail).day:
		    fig.suptitle(pd.to_datetime(str(time_head)).strftime('%B %d %Y'), fontsize=15, y=0.9)
		else:
		    fig.suptitle(pd.to_datetime(time_head).strftime('%B %d %Y') + '--' + pd.to_datetime(time_tail).strftime('%B %d %Y'),y=0.9, fontsize=15)

		choice = plt.cm.jet #seismic
		cbar_ax1 = fig.add_axes([0.91,0.305,0.01,0.095]) #[0.91,0.374,0.01,0.11]
		cbar_ax2 = fig.add_axes([0.91,0.207,0.01,0.095]) #[0.91,0.242,0.01,0.11]
		cbar_ax3 = fig.add_axes([0.91,0.11,0.01,0.095]) #[0.91,0.11,0.01,0.11]
		levels = np.linspace(-1,1,21)


		# MAGNETIC HELICITY
		ax6 = host_subplot(816)
		cs1 = ax6.contourf(Time[idx1:idx2], scale/60, sig_m, cmap=choice, levels=levels)
		cbar1 = fig.colorbar(cs1, cax=cbar_ax6)
		cbar1.ax.set_ylabel('$\\sigma_m$', fontsize=15)
		ax6.contour(Time[idx1:idx2], scale/60, sig_m, [-0.75,0.75], colors='white')
		ax6.set_ylabel('(f) \n Scale [min]', fontsize=15)
		ax6.plot(Time[idx1:idx2],coi/60,'black')
		ax6.xaxis.set_major_formatter(formatter)
		ax6.set_yscale('log',base=2)
		ax6.set_ylim(bottom=1,top=np.floor(np.max(scale)/60))
		ax6.set_xlim([time_head,time_tail])
		ax6.axes.get_xaxis().set_visible(False)
		ax6.yaxis.set_major_formatter(ticker.ScalarFormatter())
		ax6.invert_yaxis()

		# CROSS HELICITY
		ax7 = host_subplot(817)
		cs2 = ax7.contourf(Time[idx1:idx2], scale/60, sig_c, cmap=choice,levels=levels)
		cbar2 = fig.colorbar(cs2, cax=cbar_ax7)
		cbar2.ax.set_ylabel('$\\sigma_c$', fontsize=15)
		ax7.contour(Time[idx1:idx2], scale/60, sig_c, [-0.3,0.3], colors='white')
		ax7.set_ylabel('(g) \n Scale [min]', fontsize=15)
		ax7.plot(Time[idx1:idx2],coi/60,'black')
		ax7.xaxis.set_major_formatter(formatter)
		ax7.set_yscale('log',base=2)
		ax7.set_ylim(bottom=1,top=np.floor(np.max(scale)/60))
		ax7.set_xlim([time_head,time_tail])
		ax7.axes.get_xaxis().set_visible(False)
		ax7.yaxis.set_major_formatter(ticker.ScalarFormatter())
		ax7.invert_yaxis()

		# RESIDUAL ENERGY
		ax8 = host_subplot(818)
		cs3 = ax8.contourf(Time[idx1:idx2], scale/60, sig_r, cmap=choice, levels=levels)
		cbar3 = fig.colorbar(cs3, cax=cbar_ax8)
		cbar3.ax.set_ylabel('$\\sigma_r$', fontsize=15)
		ax8.set_ylabel('(h) \n Scale [min]', fontsize=15)
		ax8.plot(Time[idx1:idx2],coi/60,'black')
		ax8.xaxis.set_major_formatter(formatter)
		ax8.set_yscale('log',base=2)
		ax8.set_ylim(bottom=1,top=np.floor(np.max(scale)/60))
		ax8.set_xlim([time_head,time_tail])
		ax8.yaxis.set_major_formatter(ticker.ScalarFormatter())
		ax8.invert_yaxis()

		ticklabel = []
		for time_idx in time[np.linspace(idx1,idx2,num=6).astype(int)]:
			ticklabel.append(time_idx[11:16] + '\n %.2f\n %.2f\n %.2f' %(df_inre['x_gse'][time_idx], df_inre['y_gse'][time_idx], df_inre['z_gse'][time_idx]))

		ax8.set_xticks(df.index[np.linspace(idx1,idx2,num=6).astype(int)])
		ax8.set_xticklabels(ticklabel, fontsize=15)
		ax8.text(-0.05,-0.71,'HH:MM\n X ($R_E$)\n Y ($R_E$)\n Z ($R_E$)',horizontalalignment='right',transform=ax3.transAxes, fontsize=15)
		
		# MAGNETIC FIELD
		ax1 = host_subplot(811)
		p1a = ax4.plot(Time,Bx,color='black', label='$B_X$')
		p1b = ax4.plot(Time,By,color='blue', label='$B_Y$')
		p1c = ax4.plot(Time,Bz,color='red', label='$B_Z$')
		ax1.set_xlim([np.nanmin(Time),np.nanmax(Time)])
		ax1.xaxis.set_major_formatter(formatter)
		ax1.set_ylabel('(a) \n B [nT]', fontsize=15)
		ax1.axes.get_xaxis().set_visible(False)
		leg1 = ax1.legend(bbox_to_anchor=(1.06,0.5), loc='center right', handlelength=0, handletextpad=0, frameon=False, prop={'size': 12}, fontsize=15)
		colors=[p1a[0].get_color(),p1b[0].get_color(),p1c[0].get_color()]
		for color,text in zip(colors,leg4.get_texts()):
		    text.set_color(color)

		# VELOCITY
		ax2 = host_subplot(812)
		p2a = ax5.plot(Time,Vx,color='black', label='$V_X$')
		p2b = ax5.plot(Time,Vy,color='blue', label='$V_Y$')
		p2c = ax5.plot(Time,Vz,color='red', label='$V_Z$')
		ax2.set_xlim([np.nanmin(Time),np.nanmax(Time)])
		ax2.xaxis.set_major_formatter(formatter)
		ax2.set_ylabel('(b) \n V [km/s]', fontsize=15)
		ax2.axes.get_xaxis().set_visible(False)
		leg2 = ax2.legend(bbox_to_anchor=(1.06,0.5), loc='center right', handlelength=0, handletextpad=0, frameon=False, prop={'size': 12}, fontsize=15)
		colors=[p2a[0].get_color(),p2b[0].get_color(),p2c[0].get_color()]
		for color,text in zip(colors,leg5.get_texts()):
		    text.set_color(color)

		# DENSITY & TEMPERATURE
		ax3 = host_subplot(813)
		p3 = ax6.plot(Time,Np,color='black', label='Np')
		ax3b = ax6.twinx()
		#p2a= ax3.plot(Time,Tp*(1e-6),color='blue', label='T')
		p2a = ax3.plot(Time,Tp/11604.518,color='blue', label='T')
		ax3.set_xlim([np.nanmin(Time),np.nanmax(Time)])
		ax3.xaxis.set_major_formatter(formatter)
		ax3.set_ylabel('(c) \n $N_p$ [cm$^{-3}$]', fontsize=15)
		ax3a.set_ylabel('T [eV]', fontsize=15)
		ax3.axes.get_xaxis().set_visible(False)
		leg3 = ax3.legend(loc='best', handlelength=0, handletextpad=0, frameon=False, ncol=2, prop={'size': 12}, fontsize=15, labelcolor='linecolor')

		# ALFVEN SPEED MAGNITUDE & FLUCTUATIONS
		ax4 = host_subplot(814)
		ax4.plot(Time,bmag,color='black',label='$V_A$')
		ax4a = ax4.twinx()
		ax4a.plot(Time,Vmag-np.nanmean(Vmag),color='blue',label=r'$\delta V$')
		ax4.set_xlim([np.nanmin(Time),np.nanmax(Time)])
		ax4.xaxis.set_major_formatter(formatter)
		ax4.set_ylabel('(d) \n $V_A$ [km/s]', fontsize=15)
		ax4a.set_ylabel(r'|V| - <V> [km/s]', fontsize=13)
		ax4.axes.get_xaxis().set_visible(False)
		leg4 = ax4.legend(loc='best', handlelength=0, handletextpad=0, frameon=False, ncol=2, prop={'size': 12}, fontsize=15, labelcolor='linecolor')

		# PLASMA BETA
		ax5 = host_subplot(815)
		ax5.plot(Time,beta,color='black',label='$\\beta$')
		ax5.set_xlim([np.nanmin(Time),np.nanmax(Time)])
		ax5.xaxis.set_major_formatter(formatter)
		ax5.set_yscale('log',base=10)
		ax5.set_ylabel('(e) \n $\\beta$', fontsize=15)
		ax5.axes.get_xaxis().set_visible(False)

		plt.savefig(rootDir + 'Plots/timeseries_spectrograms_' + namestr  + '_MMS' + probe + '.png', bbox_inches='tight', dpi=300)
		#plt.show()
		plt.close()

	# ADVANCE WINDOW
	idx1 += 600
	idx2 = idx1 + 3*800
	idx1,idx2 = gap_checker(Time,idx1,idx2,dt)
