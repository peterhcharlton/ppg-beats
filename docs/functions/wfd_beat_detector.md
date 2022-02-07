# `WFD_BEAT_DETECTOR` - WFD PPG beat detector.
WFD_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
using the 'Wavelet Foot Delineation' beat detector

##  Inputs
+   sig : a vector of PPG values
    
+   fs  : the sampling frequency of the PPG in Hz
    
##  Outputs
+   peaks : indices of detected pulse peaks
    
+   onsets : indices of detected pulse onsets
    
##  Reference
N. J. Conn and D. A. Borkholder, 'Wavelet based photoplethysmogram foot delineation for heart rate variability applications,' in IEEE Signal Processing in Medicine and Biology Symposium. IEEE, 2013. <https://doi.org/10.1109/SPMB.2013.6736782>

##  Author
Elisa Mejía Mejía - wrote the code: 'wavelet'
Peter H. Charlton - did very little, just wrote this wrapper and made slight edits (see 'erma')

##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Version
0.1, and is still in development.

##  Licence
       This file is part of PPG-beats.
       PPG-beats is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
       PPG-beats is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
       You should have received a copy of the GNU General Public License along with PPG-beats. If not, see <https://www.gnu.org/licenses/>.
