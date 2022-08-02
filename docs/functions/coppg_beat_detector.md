# `COPPG_BEAT_DETECTOR` - COPPG PPG beat detector.
COPPG_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
using the 'Percentile Peak Detector' beat detector

##  Inputs
+   sig : a vector of PPG values
    
+   fs  : the sampling frequency of the PPG in Hz
    
##  Outputs
+   peaks : indices of detected pulse peaks
    
+   onsets : indices of detected pulse onsets
    
##  Reference
C. Orphanidou et al., 'Signal-quality indices for the electrocardiogram and photoplethysmogram: derivation and applications to wireless monitoring,' IEEE Journal of Biomedical and Health Informatics, vol. 19, no. 3, pp. 832-8, 2015. <https://doi.org/10.1109/JBHI.2014.2338351>

##  Author
Adapted by Peter Charlton from code supplied by Christina Orphanidou. With grateful thanks to Christina Orphanidou and Alexander Darrell.

##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Version
1.0

##  Source
+   Original source: <https://github.com/peterhcharlton/RRest/blob/master/RRest_v3.0/Algorithms/extract_resp_sig/feat_based_extraction/COr_peak_detector/co_ppg_peak_detector.m> , part of the RRest (<http://github.com/peterhcharlton/RRest/>) repository which is covered by the GNU public licence.
    
+   Modifications from original by PC: (i) Removed pre-processing; (ii) Hard-coded constants
    
##  License - GPL-3.0
       This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
       This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
       You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
