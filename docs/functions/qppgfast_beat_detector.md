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
Several authors have contributed to the code (see below)
Peter H. Charlton - did very little, just wrote this wrapper (which detects peaks given the onsets provided by the code)

##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Version
0.1, and is still in development.

##  Source
qppg_fast.m, from the PhysioNet Cardiovascular Signal Toolbox 1.0.0 at: <https://physionet.org/content/pcst/1.0.0/>

##  Licence
       This file is part of PPG-beats.
       PPG-beats is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
       PPG-beats is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
       You should have received a copy of the GNU General Public License along with PPG-beats. If not, see <https://www.gnu.org/licenses/>.