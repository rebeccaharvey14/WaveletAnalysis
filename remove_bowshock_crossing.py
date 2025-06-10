import sys
import pandas as pd
###########################################################################
rootDir = '/home/rharvey/Documents/Research/Wavelet-Analysis'

namestr = sys.argv[1]
probe = sys.argv[2]
bs_start = sys.argv[3]
bs_end = sys.argv[4]
year = sys.argv[5]

if probe.isalpha():
    filename = rootDir + 'themis_events/' + namestr[:-1] + '_THM' + probe.upper() + '.csv'
else:
    filename = rootDir + 'mms_events/' + namestr[:-1] + '_MMS' + probe + '.csv'

# EXCLUDE BOW SHOCK CROSSING
eventFile = pd.read_csv(filename,index_col=0)
print('Excluding the time between a bow shock crossing...')
before_DF = eventFile.where(eventFile['end'] < bs_start).dropna(subset=['start', 'end', 'duration', 'B_avg', 'B_max', 'beta_avg', 'V_avg', 'T_avg', 'Np_avg', 'scale_length', 'peak_time', 'sigm', 'sigc', 'sigr'])
after_DF = eventFile.where(eventFile['start'] > bs_end).dropna(subset=['start', 'end', 'duration', 'B_avg', 'B_max', 'beta_avg', 'V_avg', 'T_avg', 'Np_avg', 'scale_length', 'peak_time', 'sigm', 'sigc', 'sigr'])

print(f'Saving new wavelet event list with bow shock crossings removed to {filename}...')
eventFile = pd.concat([before_DF,after_DF],ignore_index=True)
eventFile.to_csv(filename)