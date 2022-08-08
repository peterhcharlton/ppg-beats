# `DETECT_ECG_BEATS` - detects beats in ECG.
DETECT_ECG_BEATS detects beats in an electrocardiogram (ECG) signal
using two beat detectors, and uses the agreement between beat detectors to
assess the quality of beat detections.

##  Inputs
+   ecg : a vector of ECG values
    
+   fs : the sampling frequency of the ECG in Hz
    
+   options : (optional) a structure of options which determine the settings used for the analysis:
    
     - options.win_durn       - The duration of windows over which to perform quality assessment (in secs)
     - options.win_overlap    - The overlap of windows (in secs)
     - options.verbose        - A logical indicating whether (1) or not (0) to display text in the Command Window
     - options.qrs_tol_window - The acceptable tolerance window around QRS spikes (in secs)
     - options.beat_detectors - The ECG beat detectors to be used. Provide a cell containing either 1 or 2 out of: {'gqrs', 'jqrs', 'rpeakdetect'}
     - options.hpf            - A logical indicating whether (1) or not (0) to high-pass filter the ECG to remove baseline wander
    
+   no_signal_vector : (optional) a vector of the same length as the ECG vector, which indicates whether (1) or not (0) there was a signal at each sample
    
##  Outputs
+   agreed_beats_inds : the indices of detected beats which both beat detectors agree on
    
+   qual : a logical indicating whether (1) or not (0) beats agreed between beat detectors for each ECG sample 
    
+   temporary file - this script creates a temporary file in the current directory
    
##  Exemplary usage:
        options.win_durn = 20;
        [agreed_beats_inds, qual] = detect_ecg_beats(ecg, fs, options)
    
##  Requirements:
When using the gqrs beat detector, the Waveform Database Software Package (WFDB) for MATLAB and Octave is required, available from:
<https://doi.org/10.13026/6zcz-e163>

##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Author
Peter H. Charlton, University of Cambridge, February 2022.

##  Acknowledgment:
This script uses scripts from the PhysioNet Cardiovascular Signal Toolbox (under the GNU public licence):

+   Scripts: ExportHRVparams.m, run_qrsdet_by_seg.m, jqrs.m
    
+   Some of the parameters included in this script are taken from InitializeHRVparams.m
    
+   Link: <https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox>
    
+   Reference: Vest et al. 'An open source benchmarked toolbox for cardiovascular waveform and interval analysis'. Physiol Meas, 39, 2018. <https://doi.org/10.1088/1361-6579/aae021> 
    
It also uses the 'rpeakdetect.m' script (under the GNU public licence):

+   Link: <http://www.mit.edu/~gari/CODE/ECGtools/ecgBag/rpeakdetect.m>
    
+   Author: Gari Clifford
    
##  Version
1.0

##  License - GPL-3.0
       This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
       This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
       You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
