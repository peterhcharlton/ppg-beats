# `ASSESS_BEAT_DETECTORS` - Assess beat detectors.
ASSESS_BEAT_DETECTORS assesses the performance of multiple
photoplethysmogram (PPG) beat detectors on a dataset.

##  Inputs
+   dataset - The name of the dataset with which to assess performance, e.g.
    
           'ppg_dalia_working'
           'ppg_dalia_walking'
           'ppg_dalia_sitting'
           'capnobase'
           'mimic_af'
           'mimic_non_af'
           'wesad_baseline'
    You must have downloaded and collated the dataset in order
    to use it. Instructions on how to download and collate datasets
    are provided here: <https://peterhcharlton.github.io/info/datasets>
    
+   options - (optional) A structure of options which determine the settings used for the analysis:
    
           options.dataset_file    - a string containing the path of the Matlab file containing the collated dataset.
           options.beat_detectors  - a cell array containing strings corresponding to the names of the beat detectors to be used.
           options.do_downsample   - (1 or 0, default value of 0) A logical indicating whether or not to downsample the PPG signals prior to analysis.
           options.downsample_freq - the frequency at which to downsample PPG signals (Hz). Only used if downsampling is enabled by options.do_downsample.
           options.redo_analysis   - A logical indicating whether or not to overwrite existing analysis files in order to redo the analysis
           options.redo_selected_beat_detectors - A cell containing the beat detectors for which to redo the analysis
        
##  Outputs
+   ...
    
##  Exemplary usage
        assess_beat_detectors('mimic_non_af')      % assesses the performance of beat detectors on the 'mimic_non_af' dataset using default analysis options
    
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
