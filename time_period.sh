#!/bin/bash

time_head='2009-06-19 22:00'
time_tail='2009-06-21 00:00'
namestr='_20090619_20090621'
bs_start='2009-06-20 04:45'
bs_end='2009-06-20 10:00'
probe_str='_THMC'

####################### Download And Preprocess Data #######################
############################################################################
python download_preprocess.py $time_head $time_tail $namestr $probe_str


####################### Spectrograms & Time Series Plots ###################
############################################################################
python timeseries.py $time_head $time_tail $namestr $probe_str
python timeseries_bowshock.py $time_head $time_tail $namestr $probe_str $bs_start $bs_end
python spectrogram_hourly.py $time_head $time_tail $namestr $probe_str


####################### Event Identification ###############################
############################################################################
python events.py $time_head $time_tail $namestr $probe_str


################ Remove Bow Shock Crossing (if Applicable) #################
############################################################################
python remove_bowshock_crossing.py $namestr $probe_str $bs_start $bs_end


# Save detection output files
origin_dir=/home/rharvey/Documents/Research/Wavelet-Analysis/Events/wavelet_events$namestr$probe_str.csv
destination_dir=/home/rharvey/Events/wavelet_events$namestr$probe_str.csv
echo Saving final event list file to $destination_dir
cp $origin_dir $destination_dir