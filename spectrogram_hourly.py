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
from functions import nan_checker, gap_checker, get_variables
###########################################################################

time_range = [sys.argv[1] + ' ' + sys.argv[2], sys.argv[3] + ' ' + sys.argv[4]]
namestr    = sys.argv[5]
probe_str  = sys.argv[6]
probe      = probe_str[4].lower()

rootDir = '/home/rharvey/Documents/Wavelet-Analysis/'
dataFile = '/home/rharvey/data/' + 'data' + probe_str + namestr + '.csv'

# Read in data file
Bx, By, Bz, Vx, Vy, Vz, Np, Bmag, bmag, Vmag, Tp, Te, beta, Time, dt = get_variables(dataFile,time_range)
df = pd.read_csv(dataFile,index_col=0)
df_inre = df[['x_gse','y_gse','z_gse']]
time = np.array([t.strftime('%Y-%m-%d %H:%M:%S.%f') for t in Time])
Br, Bt, Bn = -Bx, -By, Bz
Vr, Vt, Vn = -Vx, -Vy, Vz

idx1 = 0
idx2 = idx1 + 2400
while idx1 < Time.size-2400 and idx2 < Time.size:

	idx1,idx2 = gap_checker(Time,idx1,idx2,dt)
	time_head = Time[idx1]
	time_tail = Time[idx2]

	# CHECK FOR GAPS
	if nan_checker(Br[idx1:idx2]) < 0.05 and nan_checker(Vr[idx1:idx2]) < 0.05 and nan_checker(Np[idx1:idx2]) < 0.05:

		# CALCULATE MAGNETIC HELICITY, CROSS HELICITY, AND RESIDUAL ENERGY
		scale, coi, sig_m, sig_c, sig_r = helicity(Br[idx1:idx2], Bt[idx1:idx2], Bn[idx1:idx2], Vr[idx1:idx2], Vt[idx1:idx2], Vn[idx1:idx2], Np[idx1:idx2], dt)

		# PLOT OF B-FIELD AND NORMALIZED MAGNETIC HELICITY, CROSS HELICITY, & RESIDUAL ENERGY
		fig = plt.figure(figsize=(10,14))
		fig.subplots_adjust(hspace=0)
		formatter = matplotlib.dates.DateFormatter('%H:%M')
		if pd.to_datetime(time_head).day == pd.to_datetime(time_tail).day:
		    fig.suptitle(pd.to_datetime(str(time_head)).strftime('%B %d %Y'), fontsize=15, y=0.9)
		else:
		    fig.suptitle(pd.to_datetime(time_head).strftime('%B %d %Y') + '--' + pd.to_datetime(time_tail).strftime('%B %d %Y'),y=0.9, fontsize=15)

		# MAGNETIC FIELD
		ax1 = host_subplot(811)
		p1a = ax1.plot(Time[idx1:idx2],Bx[idx1:idx2],color='black', label='$B_X$')
		p1b = ax1.plot(Time[idx1:idx2],By[idx1:idx2],color='blue', label='$B_Y$')
		p1c = ax1.plot(Time[idx1:idx2],Bz[idx1:idx2],color='red', label='$B_Z$')
		ax1.set_xlim([time_head,time_tail])
		ax1.xaxis.set_major_formatter(formatter)
		ax1.set_ylabel('(a) \n B [nT]', fontsize=15)
		ax1.axes.get_xaxis().set_visible(False)
		leg1 = ax1.legend(bbox_to_anchor=(1.06,0.5), loc='center right', handlelength=0, handletextpad=0, frameon=False, prop={'size': 12}, fontsize=15, labelcolor='linecolor')

		# VELOCITY
		ax2 = host_subplot(812)
		p2a = ax2.plot(Time[idx1:idx2],Vx[idx1:idx2],color='black', label='$V_X$')
		p2b = ax2.plot(Time[idx1:idx2],Vy[idx1:idx2],color='blue', label='$V_Y$')
		p2c = ax2.plot(Time[idx1:idx2],Vz[idx1:idx2],color='red', label='$V_Z$')
		ax2.set_xlim([time_head,time_tail])
		ax2.xaxis.set_major_formatter(formatter)
		ax2.set_ylabel('(b) \n V [km/s]', fontsize=15)
		ax2.axes.get_xaxis().set_visible(False)
		leg2 = ax2.legend(bbox_to_anchor=(1.06,0.5), loc='center right', handlelength=0, handletextpad=0, frameon=False, prop={'size': 12}, fontsize=15, labelcolor='linecolor')

		# DENSITY & TEMPERATURE
		ax3 = host_subplot(813)
		p3 = ax3.plot(Time[idx1:idx2],Np[idx1:idx2],color='black', label='Np')
		ax3a = ax3.twinx()
		p2a = ax3.plot(Time[idx1:idx2],Tp[idx1:idx2]/11604.518,color='blue', label='T')
		ax3.set_xlim([time_head,time_tail])
		ax3.xaxis.set_major_formatter(formatter)
		ax3.set_ylabel('(c) \n $N_p$ [cm$^{-3}$]', fontsize=15)
		ax3a.set_ylabel('T [eV]', fontsize=15)
		ax3.axes.get_xaxis().set_visible(False)
		leg3 = ax3.legend(loc='best', handlelength=0, handletextpad=0, frameon=False, ncol=2, prop={'size': 12}, fontsize=15, labelcolor='linecolor')

		# ALFVEN SPEED MAGNITUDE & FLUCTUATIONS
		ax4 = host_subplot(814)
		ax4.plot(Time[idx1:idx2],bmag[idx1:idx2],color='black',label='$V_A$')
		ax4a = ax4.twinx()
		ax4a.plot(Time[idx1:idx2],Vmag[idx1:idx2]-np.nanmean(Vmag[idx1:idx2]),color='blue',label=r'$\delta V$')
		ax4.set_xlim([time_head,time_tail])
		ax4.xaxis.set_major_formatter(formatter)
		ax4.set_ylabel('(d) \n $V_A$ [km/s]', fontsize=15)
		ax4a.set_ylabel(r'|V| - <V> [km/s]', fontsize=13)
		ax4.axes.get_xaxis().set_visible(False)
		leg4 = ax4.legend(loc='best', handlelength=0, handletextpad=0, frameon=False, ncol=2, prop={'size': 12}, fontsize=15, labelcolor='linecolor')

		# PLASMA BETA
		ax5 = host_subplot(815)
		ax5.plot(Time[idx1:idx2],beta[idx1:idx2],color='black',label='$\\beta$')
		ax5.set_xlim([time_head,time_tail])
		ax5.xaxis.set_major_formatter(formatter)
		ax5.set_yscale('log',base=10)
		ax5.set_ylabel('(e) \n $\\beta$', fontsize=15)
		ax5.axes.get_xaxis().set_visible(False)

		choice = plt.cm.jet
		cbar_ax1 = fig.add_axes([0.91,0.305,0.01,0.095])
		cbar_ax2 = fig.add_axes([0.91,0.207,0.01,0.095])
		cbar_ax3 = fig.add_axes([0.91,0.11,0.01,0.095])
		levels = np.linspace(-1,1,21)
    
		# MAGNETIC HELICITY
		ax6 = host_subplot(816)
		cs1 = ax6.contourf(Time[idx1:idx2], scale/60, sig_m, cmap=choice, levels=levels)
		cbar1 = fig.colorbar(cs1, cax=cbar_ax1)
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
		cbar2 = fig.colorbar(cs2, cax=cbar_ax2)
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
		cbar3 = fig.colorbar(cs3, cax=cbar_ax3)
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
		ax8.text(-0.05,-0.71,'HH:MM\n X ($R_E$)\n Y ($R_E$)\n Z ($R_E$)',horizontalalignment='right',transform=ax8.transAxes, fontsize=15)

		plt.savefig(rootDir + 'Plots/timeseries_spectrograms_' + pd.to_datetime(time_head).strftime('%Y%m%d_%H%M') + probe_str + '.png', bbox_inches='tight', dpi=300)
		#plt.show()
		plt.close()

	# ADVANCE WINDOW
	idx1 += 600
	idx2 = idx1 + 2400
	idx1,idx2 = gap_checker(Time,idx1,idx2,dt)