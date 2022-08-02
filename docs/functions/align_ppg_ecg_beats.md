# `ALIGN_PPG_ECG_BEATS` - aligns PPG and ECG beats.
ALIGN_PPG_ECG_BEATS aligns beats detected in an electrocardiogram (ECG)
signal with those detected in a simultaneous photoplethysmogram (PPG) signal.

##  Inputs
+   ppg_beats_inds  - indices of PPG beats
    
+   ecg_beats_inds  - indices of ECG beats
    
+   ppg_fs          - sampling frequency of PPG in Hz
    
+   ecg_fs          - sampling frequency of ECG in Hz
    
+   options         - (optional) A structure containing the settings for the analysis:
    
     - options.max_lag     - max permissible lag between ECG and PPG in secs
     - options.lag_int     - increment between tested lags
     - options.tol_window  - acceptable tolerance between ECG and PPG beats
     - options.do_wins     - whether or not to use windowing to find the time lag
+   ecg_exc_log     - (optional) A logical indicating whether or not each ECG sample will be excluded from the analysis
    
##  Outputs
+   ecg_beats_a_inds : the aligned indices of ECG beats
    
+   ecg_exc_a_log :    the aligned ECG exclusion logical
    
+   lag_ecg_samps :    the number of samples by which the ECG is in front / behind the PPG
    
##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Author
Peter H. Charlton, University of Cambridge, February 2022.

##  Version
1.0

##  License - GPL-3.0
       Copyright (c) 2022 Peter H. Charlton
       This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
       This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
       You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
