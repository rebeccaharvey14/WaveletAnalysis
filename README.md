# WaveletAnalysis
Wavelet analysis as a method to identify structures with enhanced magnetic helicity

## File descriptions
1. **data_*.py**: downloads & preprocesses data using PySPEDAS
2. **events_*.py**: searches for events with high magnetic helicity via wavelet analysis
3. **timeseries_*.py**: plots the time series data for one observation inverval
4. **spectrogram*.py**: plots time series data + wavelet spectrograms for MHD quantities
⋅⋅* Note: **spectrogram*_hourly.py**: plots in 3-hour windows
5. **remove_bowshock_crossing.py**: removes time periods from event lists that are during a bow shock crossing (can be used for non-bowshock crossing data removal as well)
6. **MMS_time_periods.sh** and **THM_time_periods.sh**: shell scripts to implement all of the routines above
7. **helicity_calculations.py**: calculations of normalized magnetic helicity, cross helicity, and residual energy using wavelet transforms
8. **functions.py**: supplementary functions for evaluating data validity and data parsing
9. *2012 Telloni - Wavelet Analysis as a Tool to Localize Magnetic and Cross-Helicity Events in the Solar Wind.pdf* and *2021 Zhao - Detection of small magnetic flux ropes from the third and fourth Parker Solar Probe encounters.pdf*: Papers relevant to the theory and development of the algorithm implemented