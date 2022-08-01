# MIMIC III 'AF' Subset

---

This subset contains ECG and PPG recordings of 20-minute duration, some of which were acquired during atrial fibrillation (AF), and the rest were acquired during normal sinus rhythm.

## Overview

 Item | Details 
 :--- | :--- 
 **Links** | _TBC_
 **Signals** | PPG, ECG, respiration
 **No. Subjs** | 35 (19 in AF, 16 not in AF)
 **Protocol** | 35 critically-ill adults during routine clinical care. Data were measured using a bedside monitor at 125 Hz. Manual labels of AF and non-AF subjects were obtained from [here](https://doi.org/10.6084/m9.figshare.12149091.v1), as described in [[1](https://doi.org/10.1109/ACCESS.2019.2926199)]. Data were extracted from the MIMIC-III Waveform Database [Matched Subset](https://physionet.org/content/mimic3wdb-matched/1.0/).

## Importing the data into MATLAB

1. Download the MATLAB script for importing the data from [here](https://github.com/peterhcharlton/info/blob/master/collating_datasets/download_and_collate_mimiciii_ppg_af_dataset.m).
2. Modify the MATLAB script by inserting the `up.paths.root_folder` and `up.paths.save_folder` into the `setup_up` function.
   - _Note:_ The duration of recordings provided by the script can be adjusted with the `up.settings.req_durn` variable within the `setup_up` function.
3. Run the MATLAB script to download the required files from PhysioNet, import the data into MATLAB, and collate the data into a single MATLAB data file, ready for analysis.
 
