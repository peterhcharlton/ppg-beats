# `QPPGFAST_BEAT_DETECTOR` - QPPGFAST PPG beat detector.
QPPGFAST_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
using the fast version of the 'Adapted Onset Detector' beat detector

##  Inputs
+   sig : a vector of PPG values
    
+   fs  : the sampling frequency of the PPG in Hz
    
##  Outputs
+   peaks : indices of detected pulse peaks
    
+   onsets : indices of detected pulse onsets
    
##  Reference
A. N. Vest et al., 'An open source benchmarked toolbox for cardiovascular waveform and interval analysis,' Physiological Measurement, vol. 39, no. 10, 2018. <https://doi.org/10.1088/1361-6579/aae021>

##  Author
+   Several authors have contributed to the code (see below)
    
+   Peter H. Charlton - did very little, just wrote this wrapper (which detects peaks given the onsets provided by the code)
    
##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Version
1.0

##  Source
qppg_fast.m, from the PhysioNet Cardiovascular Signal Toolbox (as of 4 May 2022) at: 
<https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox>. 
The toolbox is licensed under the BSD-3-Clause License.

##  Licence
       Please see BSD 3-Clause License and GNU GPL Licence below
