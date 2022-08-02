# `SWT_BEAT_DETECTOR` - SWT PPG beat detector.
SWT_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
using the 'Stationary Wavelet Transform' beat detector

##  Inputs
+   sig : a vector of PPG values
    
+   fs  : the sampling frequency of the PPG in Hz
    
##  Outputs
+   onsets : indices of detected pulse onets
    
+   peaks : indices of detected pulse peaks
    
##  Reference (algorithm)
S. Vadrevu and M.﻿Sabarimalai Manikandan, 'A robust pulse onset and peak detection method for automated PPG signal analysis system,' IEEE Trans Instrum Meas, vol.68, no.3, pp.807-817, 2019. <https://doi.org/10.1109/TIM.2018.2857878>

##  Reference (implementation)
D. Han et al., 'A Real-Time PPG Peak Detection Method for Accurate Determination of Heart Rate during Sinus Rhythm and Cardiac Arrhythmia,' Biosensors, vol.12, no.2, p.82, 2022. <https://doi.org/10.3390/bios12020082>

##  Author
+   Dong Han - wrote the 'my_Vadrevu_2019_peakdet' function.
    
+   Peter H. Charlton - did very little, just wrote this wrapper and made a few changes to how the zero-padding of the signal is performed (as commented in the code).
    
##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Version
1.0

##  Source
my_Vadrevu_2019_peakdet.m, from the PPG_Peak_Detection GitHub Repository (accessed on 5 Mar 2022) at: <https://github.com/Cassey2016/PPG_Peak_Detection>

##  Licence
       MIT Licence (see the licence at the top of the 'my_Vadrevu_2019_peakdet' function in the code below).
