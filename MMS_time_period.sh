#!/bin/bash

probe='1'
year='2015'
time_head='2017-11-05 03:30'
time_tail='2017-11-05 12:00'
namestr='_20171105_20171106'
# bs_start=
# bs_end=

######################## Download & Preprocess Data ########################
############################################################################
#python data_MMS.py $time_head $time_tail $namestr $probe $year


####################### Spectrograms & Time Series Plots ###################
############################################################################
# python timeseries_MMS.py $time_head $time_tail $namestr $probe $year
# python spectrogram_MMS_hourly.py $time_head $time_tail $namestr $probe $year


####################### Event Identification ###############################
############################################################################
python events_MMS.py $namestr $probe


################ Remove Bow Shock Crossing (if Applicable) #################
############################################################################
#python remove_bowshock_crossing.py $namestr $probe $bs_start $bs_end $year