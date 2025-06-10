#!/bin/bash

#$1 = time_head
#$2 = time_tail
#$3 = namestr
#$4 = probe

# WAVELET SPECTROGRAMS + TIME SERIES
python3 spectrogram_MMS_hourly.py $1 $2 $3 $4

################## EVENT IDENTIFICATION ###############################
# GET LIST OF EVENTS
python3 events_MMS_v2.py $1 $2 $3 $4
