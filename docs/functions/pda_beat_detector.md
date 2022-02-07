# `PDA_BEAT_DETECTOR` - PDA PPG beat detector.
PDA_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
using the 'Peak Detection Algorithm' beat detector

##  Inputs
+   sig : a vector of PPG values
    
+   fs  : the sampling frequency of the PPG in Hz
    
##  Outputs
+   peaks : indices of detected pulse peaks
    
+   onsets : indices of detected pulse onsets
    
##  Reference
E. J. Arguello Prada and R. D. Serna Maldonado, 'A novel and low-complexity peak detection algorithm for heart rate estimation from low-amplitude photoplethysmographic (PPG) signals,' Journal of Medical Engineering and Technology, vol. 42, no. 8, pp. 569-577, 2018. <https://doi.org/10.1080/03091902.2019.1572237>

##  Author
Elisa Mejía Mejía - wrote the code: 'upslopes'
Peter H. Charlton - did very little, just wrote this wrapper

##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Version
0.1, and is still in development.

##  Licence
       This file is part of PPG-beats.
       PPG-beats is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
       PPG-beats is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
       You should have received a copy of the GNU General Public License along with PPG-beats. If not, see <https://www.gnu.org/licenses/>.
