#!/bin/bash

time_head='2009-06-19 22:00'
time_tail='2009-06-21 00:00'
namestr='_20090619_20090621'
bs_start='2009-06-20 04:45'
bs_end='2009-06-20 10:00'
probe='c'
year='2009'

####################### Download And Preprocess Data #######################
############################################################################
python data_THM.py $time_head $time_tail $namestr $probe


####################### Spectrograms & Time Series Plots ###################
############################################################################
python timeseries_THM.py $time_head $time_tail $namestr $probe
python timeseries_THM_bowshock.py $time_head $time_tail $namestr $probe $bs_start $bs_end
python spectrogram_THM_hourly.py $time_head $time_tail $namestr $probe


####################### Event Identification ###############################
############################################################################
python events_THM.py $time_head $time_tail $namestr $probe


################ Remove Bow Shock Crossing (if Applicable) #################
############################################################################
python remove_bowshock_crossing.py $namestr '_THM'$probe $bs_start $bs_end


# Save detection output files
origin_dir=/home/rharvey/Documents/Wavelet-Analysis/Events/wavelet_events$namestr_THM$probe.csv
destination_dir=/home/rharvey/Events/wavelet_events$namestr_THM$probe.csv
echo Saving final event list file to $destination_dir
cp $origin_dir $destination_dir