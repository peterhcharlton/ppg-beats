# `ASSESS_MULTIPLE_DATASETS` - Assess on multiple datasets.
ASSESS_MULTIPLE_DATASETS assesses the performance of multiple
PPG beat detectors across multiple datasets.

##  Inputs
+   none
    
+   although the datasets to be analysed should be specified in the 'setup_universal_params' function, and their paths should be specified in the 'specify_path_of_dataset_file' function.
    
##  Outputs
+   Plots illustrating the performance of PPG beat detectors.
    
+   NB: the 'ppg_beat_detector_assessment' function, which is called by this function, creates its own outputs too
    
##  Exemplary usage
        run_ppg_beat_detector_assessment      % runs the PPG beat detector assessment
    
##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Author
Peter H. Charlton, University of Cambridge, August 2022.

##  Version
1.0

##  License - GPL-3.0
       Copyright (c) 2022 Peter H. Charlton
       This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
       This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
       You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
