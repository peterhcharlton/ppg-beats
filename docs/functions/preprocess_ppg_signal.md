# `PREPROCESS_PPG_SIGNAL` - pre-processes PPG signal.
PREPROCESS_PPG_SIGNAL filters and downsamples a
photoplethysmogram (PPG) signal ready for beat detection.

##  Inputs
+   S : a structure containing three fields:
    
     - v : a vector of PPG values
     - fs : the sampling frequency of the PPG in Hz
     - no_signal : (optional) a logical indicating whether (1) or not (0) each PPG sample contains a physiological signal
    
+   options - A structure containing the settings used for the processing:
    
     - options.do_downsample   - (1 or 0, default value of 0) A logical indicating whether or not to downsample the PPG signals prior to analysis.
     - options.downsample_freq - the frequency at which to downsample PPG signals (Hz). Only used if downsampling is enabled by options.do_downsample.
     - options.redo_analysis   - A logical indicating whether or not to overwrite existing analysis files in order to redo the analysis
     - options.filtering.elim_high_freqs - ...TBC...
     - options.filtering.elim_low_freqs  - ...TBC...
    
##  Outputs
+   S : containing the filtered and downsampled PPG signal
    
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
