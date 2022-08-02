# `MIMIC_PERFORM_DATASET_CONVERTER` - converts datasets into CSV and WFDB formats.
MIMIC_PERFORM_DATASET_CONVERTER converts MIMIC PERform datasets from
Matlab format into CSV and WFDB formats.

##  Inputs
+   specify the path of the dataset for conversion in the "up.matlab_dataset.path" variable within "setup_params" below.
    
##  Outputs
+   files : containing all the data written to new subfolders within the path specified by "up.matlab_dataset.path".
    
##  Requirements:
The WFDB Toolbox is required to convert to WFDB format, which can be
downloaded from:
<https://physionet.org/physiotools/matlab/wfdb-app-matlab/>

##  Further Information
This version is provided to convert the MIMIC PERform datasets originally
described in:
Charlton P.H. et al. Detecting beats in the photoplethysmogram: benchmarking open-source algorithms, Physiological Measurement, 2022. <https://doi.org/10.1088/1361-6579/ac826d>

Further information on this study can be obtained at:
<https://peterhcharlton.github.io/project/ppg-beats/>

In addition, further information on the ppg-beats toolbox, including 
future versions, can be obtained at:
<https://github.com/peterhcharlton/ppg-beats>

##  Feedback
Questions, comments, criticisms, complements, and contributions are welcome at:
<https://github.com/peterhcharlton/ppg-beats/issues>

##  Version
v.1 - published on 1st August 2022 by Peter Charlton

##  Source
Adapted from: RRest_dataset_converter.m, part of the RRest (<https://github.com/peterhcharlton/RRest/>) repository which is covered by the GNU public licence.

##  Author
Peter H. Charlton, University of Cambridge, June 2022.

##  License - GPL-3.0
       Copyright (c) 2022 Peter H. Charlton
       This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
       This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
       You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
