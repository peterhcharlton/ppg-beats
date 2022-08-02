# `ATMIN_BEAT_DETECTOR` - ATmin PPG beat detector.
ATMIN_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
using the 'Adaptive Threshold Method (Vmin)' beat detector

##  Inputs
+   sig : a vector of PPG values
    
+   fs  : the sampling frequency of the PPG in Hz
    
##  Outputs
+   peaks : indices of detected pulse peaks
    
+   onsets : indices of detected pulse onsets
    
##  Reference (algorithm)
H. S.﻿Shin et al., 'Adaptive threshold method for the peak detection of photoplethysmographic waveform,' Comput Biol Med Elsevier, vol.39, no.12, pp.1145-52, 2009. <﻿https://doi.org/10.1016/j.compbiomed.2009.10.006>

##  Reference (implementation)
D. Han et al., 'A Real-Time PPG Peak Detection Method for Accurate Determination of Heart Rate during Sinus Rhythm and Cardiac Arrhythmia,' Biosensors, vol.12, no.2, p.82, 2022. <﻿https://doi.org/10.3390/bios12020082>

##  Author
+   Dong Han - wrote the 'my_peak_compare_Shin_2009' function.
    
+   Peter H. Charlton - did very little, just wrote this wrapper (which detects peaks given the onsets provided by the code)
    
##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Version
1.0

##  Source
my_peak_compare_Shin_2009.m, from the PPG_Peak_Detection GitHub Repository (accessed on 3 Mar 2022) at: <https://github.com/Cassey2016/PPG_Peak_Detection>

##  Licence
       MIT Licence (see the licence at the top of the 'my_peak_compare_shin_2009' function in the code below).
