# `IMS_BEAT_DETECTOR` - IMS PPG beat detector.
IMS_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
using the 'Incremental Merge Segmentation' beat detector

##  Inputs
+   sig : a vector of PPG values
    
+   fs  : the sampling frequency of the PPG in Hz
    
##  Outputs
+   peaks : detected pulse peaks
    
+   onsets : detected pulse onsets
    
##  Reference
W. Karlen et al., 'Adaptive pulse segmentation and artifact detection in photoplethysmography for mobile applications,' in Proc. IEEE EMBS. IEEE, 2012, pp. 3131-4. <https://doi.org/10.1109/EMBC.2012.6346628>

##  Author
+   Marco A. Pimentel - wrote the code: 'adaptPulseSegment'
    
+   Peter H. Charlton - did very little, just wrote this wrapper
    
##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Version
1.0

##  Source
Original source: adaptPulseSegment.m , part of the RRest (<http://github.com/peterhcharlton/RRest/>) repository which is covered by the GNU public licence.

##  License - GPL-3.0
       This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
       This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
       You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
