# `AMPD_BEAT_DETECTOR` - AMPD PPG beat detector.
AMPD_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
using the 'Automatic Multiscale-based Peak Detection' beat detector

##  Inputs
+   sig : a vector of PPG values
    
+   fs  : the sampling frequency of the PPG in Hz
    
+   lpf_freq : (optional) the frequency (in Hz) of low-pass filtering applied to the PPG
    
##  Outputs
+   peaks : indices of detected pulse peaks
    
+   onsets : indices of detected pulse onsets
    
##  Reference
F. Scholkmann et al., 'An efficient algorithm for automatic peak detection in noisy periodic and quasi-periodic signals,' Algorithms, vol. 5, no. 4, pp. 588-603, 2012. <https://doi.org/10.3390/a5040588>

##  Author
Peter H. Charlton, University of Cambridge, February 2022.

##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Version
0.1, and is still in development.

##  Licence
       This file is part of PPG-beats.
       PPG-beats is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
       PPG-beats is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
       You should have received a copy of the GNU General Public License along with PPG-beats. If not, see <https://www.gnu.org/licenses/>.
