import sys
import pandas as pd
###########################################################################
rootDir = '/home/rharvey/Documents/Research/Wavelet-Analysis/'

namestr    = sys.argv[1]
probe_str  = sys.argv[2]
bs_start   = sys.argv[3] + ' ' + sys.argv[4]
bs_end     = sys.argv[5] + ' ' + sys.argv[6]

filename = rootDir + 'Events/wavelet_events' + namestr + probe_str + '.csv'

# Exclude Bow Shock Crossing
eventFile = pd.read_csv(filename,index_col=0)
print('Excluding the time between a bow shock crossing...')
before_DF = eventFile.where(eventFile['end'] < bs_start).dropna(subset=['start', 'end', 'duration', 'B_avg', 'B_max', 'beta_avg', 'V_avg', 'T_avg', 'Np_avg', 'scale_length', 'peak_time', 'sigm', 'sigc', 'sigr'])
after_DF = eventFile.where(eventFile['start'] > bs_end).dropna(subset=['start', 'end', 'duration', 'B_avg', 'B_max', 'beta_avg', 'V_avg', 'T_avg', 'Np_avg', 'scale_length', 'peak_time', 'sigm', 'sigc', 'sigr'])

print(f'Saving new wavelet event list with bow shock crossings removed to {filename}...')
eventFile = pd.concat([before_DF,after_DF],ignore_index=True)
eventFile.to_csv(filename)