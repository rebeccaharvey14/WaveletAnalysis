# WaveletAnalysis
Wavelet analysis as a method to identify structures with enhanced magnetic helicity

## File descriptions
1. **download_preprocess.py**: downloads and preprocesses data using PySPEDAS
2. **events.py**: searches for events with high magnetic helicity via wavelet analysis
3. **timeseries.py**: plots the time series data for one observation inverval
4. **timeseries_bowshock.py**: plots the time series data, with bow shock crossing periods highlighted, for one observation inverval
5. **spectrogram_hourly.py**: plots time series data + wavelet spectrograms for MHD quantities in 3-hour windows
6. **remove_bowshock_crossing.py**: removes time periods from event list that are during a bow shock crossing (can be used for non-bowshock crossing data removal as well)
7. **time_periods.sh**: shell script to implement all of the routines above
8. **helicity_calculations.py**: calculations of normalized magnetic helicity, cross helicity, and residual energy using wavelet transforms
9. **functions.py**: supplementary functions for evaluating data validity and data parsing
10. *2012 Telloni - Wavelet Analysis as a Tool to Localize Magnetic and Cross-Helicity Events in the Solar Wind.pdf* and *2021 Zhao - Detection of small magnetic flux ropes from the third and fourth Parker Solar Probe encounters.pdf*: Papers relevant to the theory and development of the algorithm implemented