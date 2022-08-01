# PPG-DaLiA Dataset

---

## Overview

 Item | Details 
 :--- | :--- 
 **Links** | [Dataset](https://archive.ics.uci.edu/ml/datasets/PPG-DaLiA), [Publication](https://doi.org/10.3390/s19143079) 
 **Signals** | PPG, ECG (inc. manual beat annotations), respiration, accelerometry, others 
 **No. Subjs** | 15 
 **Protocol** | Recordings during a protocol of activities of daily living (car driving, cycling, in a lunch break, sitting, stair climbing, playing table soccer, walking and working)

## Importing the data into MATLAB

1. Download from [here](https://ubicomp.eti.uni-siegen.de/home/datasets/sensors19/) or [here](https://archive.ics.uci.edu/ml/datasets/PPG-DaLiA), and unzip.
2. Download the [`convert_subject_pickle_files_to_mat.py` python script](https://raw.githubusercontent.com/peterhcharlton/resources/master/collating_datasets/convert_subject_pickle_files_to_mat.py).
3. Copy the subject's _.pkl_ file to the same folder as this python script.
4. Run the python script to convert the _.pkl_ file into a MATLAB file. (I used Spyder - a Python program - through Anaconda-Navigator.)
5. Copy this MATLAB file back to the same folder as the _.pkl_ file was obtained from.
6. Download [`collate_wesad_dataset.m` MATLAB script](https://raw.githubusercontent.com/peterhcharlton/resources/master/collating_datasets/collate_wesad_dataset.m).
7. Modify the MATLAB script by inserting the `up.paths.root_folder` and `up.paths.save_folder` into the `setup_up` function.
8. Run the MATLAB script to collate all the individual MATLAB files into a single MATLAB data file, ready for analysis.