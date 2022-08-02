# MIMIC III 'Ethnicity' Subset

---

This subset contains ECG and PPG recordings of 10-minute duration, some of which were acquired from Black adults, and some of which were acquired from White adults.

## Overview

 Item | Details 
 :--- | :--- 
 **Links** | [Publication](http://peterhcharlton.github.io/publication/assess_ppg_beat_detectors/)
 **Signals** | PPG, ECG, respiration
 **No. Subjs** | 200 (100 Black, 100 White)
 **Protocol** | 200 critically-ill adults during routine clinical care. Data were measured using a bedside monitor at 125 Hz. Data were extracted from the [MIMIC-III Waveform Database Matched Subset](https://physionet.org/files/mimic3wdb-matched/1.0/).

## Importing the data into MATLAB

1. Use the following MATLAB script to import the data: [`collate_mimic_perform_ethnicity_dataset`](/functions/collate_mimic_perform_ethnicity_dataset/).
2. Modify the MATLAB script by inserting the `up.paths.root_folder` and `up.paths.save_folder` into the `setup_up` function.
   - _Note:_ The duration of recordings provided by the script can be adjusted with the `up.settings.req_durn` variable within the `setup_up` function.
   - _Note:_ The number of subjects in each group can be adjusted with the `up.settings.no_subjs_per_ethnicity` variable within the `setup_up` function.
3. Run the MATLAB script to download the required files from PhysioNet, import the data into MATLAB, and collate the data into a single MATLAB data file, ready for analysis.
4. (optionally) convert the dataset into WFDB and/or CSV formats using this MATLAB script: [`mimic_perform_dataset_converter`](/functions/mimic_perform_dataset_converter/)
