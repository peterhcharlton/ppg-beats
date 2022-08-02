# `COLLATE_MIMIC_PERFORM_TRAIN_TEST_DATASETS` - Collates training and testing datasets.
COLLATE_MIMIC_PERFORM_TRAIN_TEST_DATASETS downloads and collates data from the MIMIC-III
Waveform Database, for PPG beat detector performance evaluation.

##  Inputs
+   none (although see 'setup_up' for the parameters to be set)
    
##  Working Files
The script downloads several files from PhysioNet, which are saved locally.

##  Outputs
+   files : MATLAB files containing the required data for PPG beat detector evaluation: 'mimic_train_a_data.mat', 'mimic_train_n_data.mat', 'mimic_train_all_data.mat', 'mimic_test_a_data.mat', 'mimic_test_n_data.mat', and 'mimic_test_all_data.mat'.
    
##  Preparation
Modify the MATLAB script by inserting the 'up.paths.root_folder' and 'up.paths.save_folder' into the setup_up function.

##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Author
Peter H. Charlton, University of Cambridge, June 2022.

##  License - GPL-3.0
       Copyright (c) 2022 Peter H. Charlton
       This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
       This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
       You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
