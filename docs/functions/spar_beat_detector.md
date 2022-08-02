# `SPAR_BEAT_DETECTOR` - SPAR PPG beat detector.
SPAR_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
using the 'Symmetric Projection Attractor Reconstruction' beat detector

##  Inputs
+   sig : a vector of PPG values
    
+   fs  : the sampling frequency of the PPG in Hz
    
##  Outputs
+   peaks : indices of detected pulse peaks
    
+   onsets : indices of detected pulse onsets
    
##  Reference
none

##  Author
+   Callum Pettit and Philip Aston - wrote the code: 'SPARcycleTimesPPG' and 'find_acl'
    
+   Peter H. Charlton - modified the code slightly (as indicated), but mainly just wrote this wrapper
    
##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Version
0.1, and is still in development.

##  Licence
       Copyright (c) 2022 Callum Pettit and Philip Aston
       ...
