#!/bin/bash

time_head='2009-06-19 22:00'
time_tail='2009-06-21 00:00'
namestr='_20090619_20090621'
bs_start='2009-06-20 04:45'
bs_end='2009-06-20 10:00'
probe='c'
year='2009'

############################################################################
####################### DOWNLOAD AND PREPROCESS DATA #######################
############################################################################
python3 data_THM.py $time_head $time_tail $namestr $probe

data_dir='/home/rharvey/data/preprocessed/'
echo Saving data file to ${data_dir}DataFrame_THM${probe^^}${namestr}.csv

############################################################################
####################### SPECTROGRAMS & TIME SERIES PLOTS ###################
############################################################################
python3 timeseries_THM.py $time_head $time_tail $namestr $probe $year
python3 timeseries_THM_bowshock.py $time_head $time_tail $namestr $probe $year
python3 spectrogram_THM_hourly.py $time_head $time_tail $namestr $probe $year

############################################################################
####################### EVENT IDENTIFICATION ###############################
############################################################################
python3 events_THM.py $namestr $probe

############################################################################
####################### REMOVE BOW SHOCK CROSSING (IF APPLICABLE) ##########
############################################################################
python3 remove_bowshock_crossing.py $namestr $probe $bs_start $bs_end $year
