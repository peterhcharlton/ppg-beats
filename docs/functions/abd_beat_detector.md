# `ABD_BEAT_DETECTOR` - ABD PPG beat detector.
ABD_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
using the 'Automatic Beat Detection' beat detector

##  Inputs
+   sig : a vector of PPG values
    
+   fs  : the sampling frequency of the PPG in Hz
    
##  Outputs
+   peaks : indices of detected pulse peaks
    
+   onsets : indices of detected pulse onsets
    
##  Reference
Aboy M et al., An automatic beat detection algorithm for pressure signals. IEEE Trans Biomed Eng 2005; 52: 1662-70. <https://doi.org/10.1109/TBME.2005.855725>

##  Author
Peter H. Charlton: King's College London (August 2017), University of Cambridge (February 2022)

##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Version
1.0

##  Source
This script contains items either copied or modified from the pulse-analyse
(<https://github.com/peterhcharlton/pulse-analyse>) and RRest 
(<http://github.com/peterhcharlton/RRest/>) repositories which are covered
by the GNU public licence.

##  License - GPL-3.0
       Copyright (c) 2022 Peter H. Charlton
       This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
       This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
       You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
