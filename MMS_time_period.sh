#!/bin/bash

probe='1'
year='2015'
time_head='2017-11-05 03:30'
time_tail='2017-11-05 12:00'
namestr='_20171105_20171106'
# bs_start=
# bs_end=

############################################################################
####################### DOWNLOAD AND PREPROCESS DATA #######################
############################################################################
python3 data_MMS.py $time_head $time_tail $namestr $probe

data_dir='/home/rharvey/data/preprocessed/'
echo Saving data file to ${data_dir}DataFrame_MMS1${namestr}.csv

############################################################################
####################### SPECTROGRAMS & TIME SERIES PLOTS ###################
############################################################################
python3 timeseries_MMS.py $time_head $time_tail $namestr $probe $year
python3 spectrogram_MMS_hourly.py $time_head $time_tail $namestr $probe $year

############################################################################
####################### EVENT IDENTIFICATION ###############################
############################################################################
python3 events_MMS.py $namestr $probe

############################################################################
####################### REMOVE BOW SHOCK CROSSING (IF APPLICABLE) ##########
############################################################################
python3 remove_bowshock_crossing.py $namestr $probe $bs_start $bs_end $year