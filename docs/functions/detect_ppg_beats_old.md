# `DETECT_PPG_BEATS` - detects beats in PPG.
DETECT_PPG_BEATS detects beats in a photoplethysmogram (PPG) signal
using a specified beat detector

##  Inputs
+   s : a structure containing the following fields:
    
     - v : a vector of PPG values
     - fs : the sampling frequency of the PPG in Hz
    
+   beat_detector  - a string specifying the beat detector to be used
    
##  Outputs
+   peaks : ...TBC...
    
+   onsets : ...TBC...
    
+   mid_amps : ...TBC...
    
##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Author
Peter H. Charlton, University of Cambridge, February 2022.

##  Version
0.1, and is still in development.

##  Licence
       This file is part of PPG-beats.
       PPG-beats is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
       PPG-beats is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
       You should have received a copy of the GNU General Public License along with PPG-beats. If not, see <https://www.gnu.org/licenses/>.
