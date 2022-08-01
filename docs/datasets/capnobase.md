# CapnoBase Dataset

---

## Overview

 Item | Details 
 :--- | :--- 
 **Links** | [Dataset](https://doi.org/10.5683/SP2/NLB8IT), [Publication](https://doi.org/10.1109/TBME.2013.2246160) 
 **Signals** | PPG, ECG (inc. manual beat annotations), respiration
 **No. Subjs** | 42 
 **Protocol** | Recordings from adults and children during elective surgery and routine anaesthesia (mostly high-quality recordings).  

## Importing the data into MATLAB

I took the following steps to import the CapnoBase dataset to MATLAB:

1. Download the data from [here](https://doi.org/10.5683/SP2/NLB8IT). Note that the dataset contains 337 files, but only the 42 files ending in __8min.mat_ are required. You can select these files by: (i) searching for "*_8min.mat" in the 'Find' box; (ii) setting the 'Files Per Page' to 50 (see the dropdown menu at the bottom of the page; and (iii) ticking the box at the header of the table to select all 42 listed files. Then use the 'Download' button on the right hand side of the header to download the files.
2. Unzip the downloaded zip file.
3. Download the [`collate_capnobase_dataset.m` MATLAB script](https://raw.githubusercontent.com/peterhcharlton/info/master/collating_datasets/collate_capnobase_dataset.m).
4. Modify the MATLAB script by inserting the `up.paths.root_folder` and `up.paths.save_folder` into the `setup_up` function.
5. Run the MATLAB script to collate all the individual MATLAB files into a single MATLAB data file, ready for analysis.
