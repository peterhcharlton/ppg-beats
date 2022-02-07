# `HEARTPY_BEAT_DETECTOR` - HEARTPY PPG beat detector.
HEARTPY_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
using the 'HeartPy' beat detector

##  Inputs
+   sig : a vector of PPG values
    
+   fs  : the sampling frequency of the PPG in Hz
    
##  Outputs
+   peaks : indices of detected pulse peaks
    
+   onsets : indices of detected pulse onsets
    
##  Reference
P. van Gent et al., 'HeartPy: A novel heart rate algorithm for the analysis of noisy signals,' Transportation Research Part F: Traffic Psychology and Behaviour, vol. 66, pp. 368-378, 2019. https://doi.org/10.1016/j.trf.2019.09.015>

##  Author
Peter H. Charlton, University of Cambridge, February 2022.

##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Version
0.1, and is still in development.

##  Source
This is a Matlab implementation of the original Python algorithm, described in:
van Gent P. et al., Analysing noisy driver physiology real-time using off-the-shelf sensors: Heart rate analysis software from the taking the fast lane project. J. Open Res. Softw. 2019, 7, <https://doi.org/10.5334/jors.241>
It was created based on v.1.2.5 of HeartPy downloaded on 7-June-2021 from: <https://github.com/paulvangentcom/heartrate_analysis_python> 
NB: The original HeartPy algorithm contains several extra features (such as artifact detection), which are not included in this code.

##  Licence
       This file is part of PPG-beats.
       PPG-beats is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
       PPG-beats is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
       You should have received a copy of the GNU General Public License along with PPG-beats. If not, see <https://www.gnu.org/licenses/>.
