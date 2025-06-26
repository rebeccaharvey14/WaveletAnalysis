import datetime
import pandas as pd
import numpy as np

def get_variables(filename,time_range):
    df = pd.read_csv(filename,index_col=0)
    time_head, time_tail = time_range
    selected_index = [(df.index >= time_head) & (df.index <= time_tail)][0]
    df = df[selected_index]

    Time = pd.to_datetime(df.index)
    dt = Time[1]-Time[0]
    dt = dt.seconds + dt.microseconds*1e-6
    Np = df.Np.values
    Tp = df.Tp.values
    Te = df.Te.values
    beta = df.beta.values
    bx = df.bx.values
    by = df.by.values
    bz = df.bz.values
    bmag = df.bmag.values
    Bx = df.Bx.values
    By = df.By.values
    Bz = df.Bz.values
    Bmag = df.Bmag.values
    Vx = df.Vx.values
    Vy = df.Vy.values
    Vz = df.Vz.values
    Vmag = df.Vmag.values
    return Bx, By, Bz, Vx, Vy, Vz, Np, Bmag, bmag, Vmag, Tp, Te, beta, Time, dt

def nan_checker(array):
	idx_nans = np.where(np.isnan(array)==True)[0]
	return idx_nans.size/array.size

def gap_checker(Time,idx1,idx2,dt):
	if idx2 > Time.size:
		idx2 = Time.size-1
	elif (Time[idx2]-Time[idx1]) > datetime.timedelta(seconds=3600*3):
		while (Time[idx2]-Time[idx1]) > datetime.timedelta(seconds=2400*dt):
			idx2 -= 1
	return idx1,idx2