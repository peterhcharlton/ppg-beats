# Datasets

This page provides details of ...

---

## CapnoBase Dataset

### Overview

 Item | Details 
 :--- | :--- 
 **Links** | [Dataset](https://doi.org/10.5683/SP2/NLB8IT), [Publication](https://doi.org/10.1109/TBME.2013.2246160) 
 **Signals** | PPG, ECG, resp
 **No. Subjs** | 42 
 **Protocol** | Recordings from adults and children during elective surgery and routine anaesthesia  

### Importing the data into MATLAB

I took the following steps to import the CapnoBase dataset to MATLAB:

1. Download the data from [here](https://doi.org/10.5683/SP2/NLB8IT). Note that the dataset contains 337 files, but only the 42 files ending in __8min.mat_ are required. You can select these files by: (i) searching for "*_8min.mat" in the 'Find' box; (ii) setting the 'Files Per Page' to 50 (see the dropdown menu at the bottom of the page; and (iii) ticking the box at the header of the table to select all 42 listed files. Then use the 'Download' button on the right hand side of the header to download the files.
2. Unzip the downloaded zip file.
3. Download the [`collate_capnobase_dataset.m` MATLAB script](https://raw.githubusercontent.com/peterhcharlton/info/master/collating_datasets/collate_capnobase_dataset.m).
4. Modify the MATLAB script by inserting the `up.paths.root_folder` and `up.paths.save_folder` into the `setup_up` function.
5. Run the MATLAB script to collate all the individual MATLAB files into a single MATLAB data file, ready for analysis.

---

## BIDMC Dataset

### Overview

 Item | Details 
 :--- | :--- 
 **Links** | [Dataset](https://doi.org/10.13026/C2208R), [Publication](https://doi.org/10.1109/TBME.2016.2613124) 
 **Signals** | PPG, ECG, resp
 **No. Subjs** | 53 
 **Protocol** | Recordings from critically-ill adults during routine clinical care. Data were measured using a bedside monitor at 125 Hz. Data were extracted from the [MIMIC II Matched Waveform Database](https://physionet.org/physiobank/database/mimic2wdb/matched/).

### Importing the data into MATLAB

1. Download the dataset in MATLAB format from [here](https://physionet.org/content/bidmc/1.0.0/bidmc_data.mat).

---

## MIMIC 'PERform' Training Dataset

This subset contains ECG and PPG recordings of 10-minute duration, some of which were acquired from adults, and some of which were acquired neonates.

### Overview

 Item | Details 
 :--- | :--- 
 **Links** | _TBC_
 **Signals** | PPG, ECG, resp
 **No. Subjs** | 200 (100 adults, 100 neonates)
 **Protocol** | 200 critically-ill patients during routine clinical care. Data were measured using a bedside monitor at 125 Hz. Data were extracted from the MIMIC-III [Waveform Database](https://physionet.org/content/mimic3wdb/1.0/).

### Importing the data into MATLAB

1. Download the MATLAB script for importing the data from _TBC_.
2. Modify the MATLAB script by inserting the `up.paths.root_folder` and `up.paths.save_folder` into the `setup_up` function.
   - _Note:_ The duration of recordings provided by the script can be adjusted with the `up.settings.req_durn` variable within the `setup_up` function.
   - _Note:_ The number of subjects in each group can be adjusted with the `up.settings.no_subjs_per_group` variable within the `setup_up` function.
3. Run the MATLAB script to download the required files from PhysioNet, import the data into MATLAB, and collate the data into a single MATLAB data file, ready for analysis.

---

## MIMIC 'PERform' Testing Dataset

This subset contains ECG and PPG recordings of 10-minute duration, some of which were acquired from adults, and some of which were acquired neonates.

### Overview

 Item | Details 
 :--- | :--- 
 **Links** | _TBC_
 **Signals** | PPG, ECG, resp
 **No. Subjs** | 200 (100 adults, 100 neonates)
 **Protocol** | 200 critically-ill patients during routine clinical care. Data were measured using a bedside monitor at 125 Hz. Data were extracted from the MIMIC-III [Waveform Database](https://physionet.org/content/mimic3wdb/1.0/).

### Importing the data into MATLAB

1. Download the MATLAB script for importing the data from _TBC_.
2. Modify the MATLAB script by inserting the `up.paths.root_folder` and `up.paths.save_folder` into the `setup_up` function.
   - _Note:_ The duration of recordings provided by the script can be adjusted with the `up.settings.req_durn` variable within the `setup_up` function.
   - _Note:_ The number of subjects in each group can be adjusted with the `up.settings.no_subjs_per_group` variable within the `setup_up` function.
3. Run the MATLAB script to download the required files from PhysioNet, import the data into MATLAB, and collate the data into a single MATLAB data file, ready for analysis.

---

## MIMIC III 'AF' Subset

This subset contains ECG and PPG recordings of 1-hour duration, some of which were acquired during atrial fibrillation (AF), and the rest were acquired during normal sinus rhythm.

### Overview

 Item | Details 
 :--- | :--- 
 **Links** | _TBC_
 **Signals** | PPG, ECG, resp
 **No. Subjs** | 34 (18 in AF, 16 not in AF)
 **Protocol** | 34 critically-ill adults during routine clinical care. Data were measured using a bedside monitor at 125 Hz. Manual labels of AF and non-AF subjects were obtained from [here](https://doi.org/10.6084/m9.figshare.12149091.v1), as described in [[1](https://doi.org/10.1109/ACCESS.2019.2926199)]. Data were extracted from the MIMIC-III Waveform Database [Matched Subset](https://physionet.org/content/mimic3wdb-matched/1.0/).

### Importing the data into MATLAB

1. Download the MATLAB script for importing the data from [here](https://github.com/peterhcharlton/info/blob/master/collating_datasets/download_and_collate_mimiciii_ppg_af_dataset.m).
2. Modify the MATLAB script by inserting the `up.paths.root_folder` and `up.paths.save_folder` into the `setup_up` function.
   - _Note:_ The duration of recordings provided by the script can be adjusted with the `up.settings.req_durn` variable within the `setup_up` function.
3. Run the MATLAB script to download the required files from PhysioNet, import the data into MATLAB, and collate the data into a single MATLAB data file, ready for analysis.
 
 ---
 
## MIMIC III 'Ethnicity' Subset

This subset contains ECG and PPG recordings of 10-minute duration, some of which were acquired from Black adults, and some of which were acquired from White adults.

### Overview

 Item | Details 
 :--- | :--- 
 **Links** | _TBC_
 **Signals** | PPG, ECG, resp
 **No. Subjs** | 100 (50 Black, 50 White)
 **Protocol** | 100 critically-ill adults during routine clinical care. Data were measured using a bedside monitor at 125 Hz. Data were extracted from the MIMIC-III [Waveform Database](https://physionet.org/content/mimic3wdb/1.0/).

### Importing the data into MATLAB

1. Download the MATLAB script for importing the data from _TBC_.
2. Modify the MATLAB script by inserting the `up.paths.root_folder` and `up.paths.save_folder` into the `setup_up` function.
   - _Note:_ The duration of recordings provided by the script can be adjusted with the `up.settings.req_durn` variable within the `setup_up` function.
   - _Note:_ The number of subjects in each group can be adjusted with the `up.settings.no_subjs_per_ethnicity` variable within the `setup_up` function.
3. Run the MATLAB script to download the required files from PhysioNet, import the data into MATLAB, and collate the data into a single MATLAB data file, ready for analysis.

---

## WESAD Dataset

### Overview

 Item | Details 
 :--- | :--- 
 **Links** | [Dataset](https://archive.ics.uci.edu/ml/datasets/WESAD+%28Wearable+Stress+and+Affect+Detection%29), [Publication](https://doi.org/10.1145/3242969.3242985) 
 **Signals** | PPG, ECG, resp, accel, others 
 **No. Subjs** | 15 
 **Protocol** | Recordings at baseline and during amusement, stress and meditation 

### Importing the data into MATLAB

I took the following steps to import the WESAD dataset to MATLAB:

1. Download the data from [here](https://uni-siegen.sciebo.de/s/HGdUkoNlW1Ub0Gx) (as stated [here](https://ubicomp.eti.uni-siegen.de/home/datasets/icmi18/)).
2. Download the [`convert_subject_pickle_files_to_mat.py` python script](https://raw.githubusercontent.com/peterhcharlton/resources/master/collating_datasets/convert_subject_pickle_files_to_mat.py).
3. Copy the subject's _.pkl_ file to the same folder as this python script.
4. Run the python script to convert the _.pkl_ file into a MATLAB file. (I used Spyder - a Python program - through Anaconda-Navigator.)
5. Copy this MATLAB file back to the same folder as the _.pkl_ file was obtained from.
6. Download the [`collate_ppg_dalia_dataset.m` MATLAB script](https://raw.githubusercontent.com/peterhcharlton/resources/master/collating_datasets/collate_ppg_dalia_dataset.m).
7. Modify the MATLAB script by inserting the `up.paths.root_folder` and `up.paths.save_folder` into the `setup_up` function.
8. Run the MATLAB script to collate all the individual MATLAB files into a single MATLAB data file, ready for analysis.

---

## PPG-DaLiA Dataset

### Overview

 Item | Details 
 :--- | :--- 
 **Links** | [Dataset](https://archive.ics.uci.edu/ml/datasets/PPG-DaLiA), [Publication](https://doi.org/10.3390/s19143079) 
 **Signals** | PPG, ECG, resp, accel, others 
 **No. Subjs** | 15 
 **Protocol** | Recordings during activities of daily living (car driving, cycling, in a lunch break, sitting, stair climbing, playing table soccer, walking and working)

### Importing the data into MATLAB

1. Download from [here](https://ubicomp.eti.uni-siegen.de/home/datasets/sensors19/) or [here](https://archive.ics.uci.edu/ml/datasets/PPG-DaLiA), and unzip.
2. Download the [`convert_subject_pickle_files_to_mat.py` python script](https://raw.githubusercontent.com/peterhcharlton/resources/master/collating_datasets/convert_subject_pickle_files_to_mat.py).
3. Copy the subject's _.pkl_ file to the same folder as this python script.
4. Run the python script to convert the _.pkl_ file into a MATLAB file. (I used Spyder - a Python program - through Anaconda-Navigator.)
5. Copy this MATLAB file back to the same folder as the _.pkl_ file was obtained from.
6. Download [`collate_wesad_dataset.m` MATLAB script](https://raw.githubusercontent.com/peterhcharlton/resources/master/collating_datasets/collate_wesad_dataset.m).
7. Modify the MATLAB script by inserting the `up.paths.root_folder` and `up.paths.save_folder` into the `setup_up` function.
8. Run the MATLAB script to collate all the individual MATLAB files into a single MATLAB data file, ready for analysis.