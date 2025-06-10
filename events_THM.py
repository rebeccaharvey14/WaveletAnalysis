# FINDS CONTOURS OF SPECTROGRAMS AND IDENTIFIES DURATION & PARAMETERS OF EVENTS

import os
import sys
import datetime
import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from waveletFunctions import wavelet
from my_functions import find_nearest
from helicity_calculations import helicity
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
####################### User defined module #######################
def get_variables(filename, rtn=True):
	df = pd.read_pickle(filename)
	dj = 0.125
	Np = df.Np.values
	Bmag = df.Bmag.values
	Tp = df.Tp.values
	Te = df.Te.values
	Time = df.index
	Vmag = df.Vmag.values
	beta = df.beta.values

	dt = Time[1]-Time[0]
	dt = dt.seconds + dt.microseconds*1e-6

	if 'rtn':
		Br = -df.Bx.values
		Bt = -df.By.values
		Bn = df.Bz.values
		Vr = -df.Vx.values
		Vt = -df.Vy.values
		Vn = df.Vz.values
		return Br, Bt, Bn, Vr, Vt, Vn, Np, Bmag, Vmag, Tp, Te, beta, Time, dj, dt
	else:
		Bx = df.Bx.values
		By = df.By.values
		Bz = df.Bz.values
		Vx = df.Vx.values
		Vy = df.Vy.values
		Vz = df.Vz.values
		return Bx, By, Bz, Vx, Vy, Vz, Np, Bmag, Vmag, Tp, Te, beta, Time, dj, dt

def nan_checker(array):
	idx_nans = np.where(np.isnan(array)==True)[0]
	return idx_nans.size/array.size

def gap_checker(array,idx1,idx2,dt):
	if idx2 == array.size or idx2>(array.size-20):
		idx2 = array.size-1
	elif (array[idx2]-array[idx1]) > datetime.timedelta(seconds=3600*3):
		while (array[idx2]-array[idx1]) > datetime.timedelta(seconds=3600*3):
			idx2 -= 1
		if idx1 == idx2 or idx2-idx1==1: # if there no points
			idx1 = idx2 + 600
			idx2 = idx1 + 2400
	return idx1,idx2
###########################################################################

# READ IN DATA
namestr = sys.argv[1]
probe = sys.argv[2]
rootDir = '/home/rharvey/Documents/Research/Wavelet-Analysis/'
datafile = '/home/rharvey/data/preprocessed/' + 'DataFrame_THM' + probe.upper() + namestr + '.csv'
eventFile = rootDir + 'themis_events/' + namestr[:-1] + '_THM' + probe.upper() + '.csv'

Br, Bt, Bn, Vr, Vt, Vn, Np, Bmag, Vmag, Tp, Te, beta, Time, dj, dt = get_variables(datafile,rtn=True)
epoch = ((Time-pd.Timestamp('1970-01-01')) // pd.Timedelta('1s')).values
Time = np.array([datetime.datetime.utcfromtimestamp(t) for t in epoch])

cols = ['start', 'end', 'duration', 'B_avg', 'B_max', 'beta_avg', 'V_avg', 'T_avg', 'Np_avg', 'scale_length', 'peak_time', 'sigm', 'sigc', 'sigr']
eventList = pd.DataFrame(columns=cols)

idx1 = 0
idx2 = idx1 + 2400
while idx1 < Time.size-2400 and idx2 < Time.size:

	time_head = Time[idx1]
	time_tail = Time[idx2]

	if nan_checker(Br[idx1:idx2]) < 0.05 and nan_checker(Vr[idx1:idx2]) < 0.05 and nan_checker(Np[idx1:idx2]) < 0.05:

		idx_nans_B = np.where(np.isnan(Br[idx1:idx2])==True)[0]
		idx_nans_V = np.where(np.isnan(Vr[idx1:idx2])==True)[0]
		idx_nans_Np = np.where(np.isnan(Np[idx1:idx2])==True)[0]

		# CALCULATE MAGNETIC HELICITY, CROSS HELICITY, AND RESIDUAL ENERGY
		scale, coi, sig_m, sig_c, sig_r = helicity(Br[idx1:idx2], Bt[idx1:idx2], Bn[idx1:idx2], Vr[idx1:idx2], Vt[idx1:idx2], Vn[idx1:idx2], Np[idx1:idx2], dt, dj,magnetosheath=True)

		# CONDITIONS: >10*dt, <3 hours
		idx_max = np.where(scale<(3*60*60))[0]
		idx_min = np.where(scale[idx_max]>10*dt)[0]
		
		if idx_min.size !=0:
			[r,c] = sig_m[idx_min,:].shape

			# FIND LOCATION OF CONTOUR INTERVALS
			fig, ax = plt.subplots(1)
			ct = ax.contour(Time[idx1:idx2], scale, sig_m, [-0.75,0.75], colors=['w'])
			arr_1 = np.asarray(ct.allsegs,dtype=object)
			
			if arr_1.shape == (2,):
			
				arr_minus = np.asarray(arr_1[0],dtype=object) # negative threshold
				arr_plus = np.asarray(arr_1[1],dtype=object) # positive threshold
				Time_num = matplotlib.dates.date2num(list(Time[idx1:idx2]))
				plt.close()

				n = len(arr_minus)
				m = len(arr_plus)
				scale_arr = []
				start = []
				end = []
				B_avg = []
				B_max = []
				V_avg = []
				beta_avg = []
				temp_avg = []
				xidx = []
				yidx = []
				max_arr = []
				Np_avg = []
				scale_length = []
				Vmag_arr = []
				peak_time = []
				coi_arr = []
				sigc = []
				sigr = []

				for i in range(n+m):
					# time, scale coordinates of contour region
					if i < n:
						[x,y] = np.hsplit(arr_minus[i],2)
					else:
						[x,y] = np.hsplit(arr_plus[i-n],2)
						
					x = x.reshape((1,x.size))[0]
					y = y.reshape((1,y.size))[0]
					sort_idx = np.argsort(x)
					x = x[sort_idx]
					y = y[sort_idx]

					idx_time = find_nearest(x,Time_num).astype(int)
					idx_scale = find_nearest(y,scale).astype(int)
					
					# IF THERE IS ANY SCALE >10*DT, <6 HOURS, BEGIN LOOP
					if any(idx_scale>idx_min[0]) and any(idx_scale<idx_min[idx_min.size-1]):

						maximum = 0
						for j,k in zip(idx_scale,idx_time):
							if j>idx_min[0] and j<idx_min[idx_min.size-1] and np.mean(scale[idx_scale])<np.mean(coi[idx_time]):
								if np.abs(sig_m[j,k]) > np.abs(maximum):
									maximum = sig_m[j,k]
									
						if maximum!=0:
							max_idx = np.where(sig_m==maximum)

							scale_arr = np.append(scale_arr, scale[max_idx[0][0]]/60) # minutes
							scale_length = np.append(scale_length, np.mean(Vmag[idx1:idx2][idx_time])*scale[max_idx[0][0]]) # km
							start = np.append(start, Time[idx1:idx2][k] - datetime.timedelta(seconds=0.5*scale[max_idx[0][0]]))
							end = np.append(end, Time[idx1:idx2][k] + datetime.timedelta(seconds=0.5*scale[max_idx[0][0]]))
							
							B_avg = np.append(B_avg, np.mean(Bmag[idx1:idx2][idx_time]))
							B_max = np.append(B_max, np.max(Bmag[idx1:idx2][idx_time]))
							V_avg = np.append(V_avg, np.mean(Vmag[idx1:idx2][idx_time]))
							beta_avg = np.append(beta_avg, np.mean(beta[idx1:idx2][idx_time]))
							temp_avg = np.append(temp_avg, np.mean(Tp[idx1:idx2][idx_time]))
							Np_avg = np.append(Np_avg, np.mean(Np[idx1:idx2][idx_time]))
							xidx = np.append(xidx,max_idx[1][0])
							yidx = np.append(yidx,max_idx[0][0])
							max_arr = np.append(max_arr, maximum)
							peak_time = np.append(peak_time, Time[idx1:idx2][k])
							#sigm = np.append(sigm, sig_m[max_idx])
							sigc = np.append(sigc, sig_c[max_idx])
							sigr = np.append(sigr, sig_r[max_idx])

				eventList_interval = pd.DataFrame(columns=cols)
				eventList_interval['sigm'] = max_arr
				eventList_interval['start']	= start 
				eventList_interval['end'] = end 
				eventList_interval['duration'] = scale_arr
				eventList_interval['B_avg'] = B_avg
				eventList_interval['B_max'] = B_max
				eventList_interval['beta_avg'] = beta_avg
				eventList_interval['V_avg'] = V_avg
				eventList_interval['T_avg'] = temp_avg
				eventList_interval['Np_avg'] = Np_avg
				eventList_interval['scale_length'] = scale_length
				eventList_interval['peak_time'] = peak_time
				eventList_interval['sigc'] = sigc
				eventList_interval['sigr'] = sigr
				eventList = pd.concat([eventList, eventList_interval],axis=0,ignore_index=True)
			plt.close()
	idx1 += 600
	idx2 = idx1 + 2400
	idx1,idx2 = gap_checker(Time,idx1,idx2,dt)

# ROUND VALUES IN THE EVENT LIST TO 3 DECIMALS
eventList = eventList.round({'sigm':3, 'duration':5, 'B_avg':5, 'B_max':5, 'beta_avg':5, 'V_avg':5, 'T_avg':5, 'Np_avg':5, 'scale_length':5, 'sigc':3, 'sigr':3})
eventList = eventList.drop_duplicates()
eventList.sort_values('start', axis=0, ascending=True, inplace=True, kind='quicksort', ignore_index=True)
# print('All events...')
# print(eventList)

# DROP EVENTS WITH THE SAME PEAK TIME
print('\nEliminating events with same peak time...')
# eventList.drop_duplicates().drop_duplicates(subset=['peak_time']).reset_index(drop=True,inplace=True)
# print(eventList)

# ELIMINATE EVENTS THAT OVERLAP
overlap = len(eventList) +1
while overlap > len(eventList):
	# CREATE ARRAY OF INTERVALS
	start_array = np.array([pd.Timestamp(elem) for elem in eventList['start'].values])
	end_array = np.array([pd.Timestamp(elem) for elem in eventList['end'].values])
	interval_array = pd.arrays.IntervalArray.from_arrays(start_array,end_array,closed='both')

	# RETRIEVING COLUMN INDEX FROM COLUMN NAME.
	eventList = eventList.assign(keepFlag=[True]*len(eventList))
	idx_column_keepFlag = eventList.columns.get_loc('keepFlag')

	# KEEP EVENT WITH HIGHEST DURATION
	for i_index in range(len(eventList)):
		interval = interval_array[i_index]
		idx_overlap = np.where(interval_array.overlaps(interval)==True)[0]
		max_duration_idx = eventList.loc[idx_overlap]['duration'].idxmax()
		eventList.iloc[idx_overlap, idx_column_keepFlag] = False
		eventList.loc[max_duration_idx, 'keepFlag'] = True

	eventList = eventList[eventList['keepFlag']==True]
	eventList = eventList.drop(columns=['keepFlag'])
	eventList.reset_index(drop=True,inplace=True)

	start_array = np.array([pd.Timestamp(elem) for elem in eventList['start'].values])
	end_array = np.array([pd.Timestamp(elem) for elem in eventList['end'].values])
	interval_array = pd.arrays.IntervalArray.from_arrays(start_array,end_array,closed='both')
	overlap = np.array([interval_array.overlaps(interval_array[i]).sum() for i in range(len(interval_array))]).sum()

noOverlap = eventList.copy()
noOverlap.sort_values('start', axis=0, ascending=True, inplace=True, kind='quicksort', ignore_index=True)

if os.path.isfile(eventFile):
	eventFile = rootDir + 'mms_events/' + namestr + '_THM' + probe.upper() + '.csv'

print(f'\nSaving event list to {eventFile}')
noOverlap.to_csv(eventFile)
print(noOverlap)
