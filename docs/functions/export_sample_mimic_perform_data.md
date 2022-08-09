# `EXPORT_SAMPLE_MIMIC_PERFORM_DATA` - Export sample data.
EXPORT_SAMPLE_MIMIC_PERFORM_DATA exports sample data excerpts from the
MIMIC Perform datasets.

##  Inputs
None, although the script uses the following MIMIC PERform dataset files (which can be downloaded from <https://doi.org/10.5281/zenodo.6807402>):

+   mimic_perform_train_all_data.mat
    
+   mimic_perform_train_n_data.mat
    
+   mimic_perform_test_all_data.mat
    
+   mimic_perform_af_data.mat
    
+   mimic_perform_non_af_data.mat
    
##  Outputs
+   files : MATLAB files containing the data excerpts.
    
##  Preparation:
Modify the MATLAB script by specifying the filepaths of the data files and the folder in which to save excerpts (within the setup_up function).

##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Author
Peter H. Charlton, University of Cambridge, June 2022.

##  License - GPL-3.0
       Copyright (c) 2022 Peter H. Charlton
       This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
       This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
       You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
