# MIMIC PERform AF Dataset

---

This dataset contains ECG and PPG recordings of 20-minute duration, some of which were acquired during atrial fibrillation (AF), and the rest were acquired during normal sinus rhythm.

## Overview

 Item | Details 
 :--- | :--- 
 **Links** | [Dataset](https://doi.org/10.5281/zenodo.6807402), [Publication](http://peterhcharlton.github.io/publication/assess_ppg_beat_detectors/)
 **Signals** | PPG, ECG, respiration
 **No. Subjs** | 35 (19 in AF, 16 not in AF)
 **Protocol** | 35 critically-ill adults during routine clinical care. Data were measured using a bedside monitor at 125 Hz. Manual labels of AF and non-AF subjects were obtained from [here](https://doi.org/10.6084/m9.figshare.12149091.v1), as described in [[1](https://doi.org/10.1109/ACCESS.2019.2926199)]. Data were extracted from the MIMIC-III Waveform Database [Matched Subset](https://physionet.org/content/mimic3wdb-matched/1.0/).
 
## Downloading the dataset

Links to download the dataset in different formats are provided in this table:
 
Format | Links 
:--- | :--- 
Matlab | [AF subjects](https://zenodo.org/record/6807403/files/mimic_perform_af_data.mat?download=1), [non-AF subjects](https://zenodo.org/record/6807403/files/mimic_perform_non_af_data.mat?download=1)
CSV | [AF subjects](https://zenodo.org/record/6807403/files/mimic_perform_af_csv.zip?download=1), [non-AF subjects](https://zenodo.org/record/6807403/files/mimic_perform_non_af_csv.zip?download=1)
WFDB | [AF subjects](https://zenodo.org/record/6807403/files/mimic_perform_af_wfdb.zip?download=1), [non-AF subjects](https://zenodo.org/record/6807403/files/mimic_perform_non_af_wfdb.zip?download=1)

## Recreating the dataset

The dataset can be recreated from scratch by taking the following steps:

1. Download the following MATLAB script: [`collate_mimic_perform_af_dataset`](https://ppg-beats.readthedocs.io/en/latest/functions/collate_mimic_perform_af_dataset/).
2. Modify the MATLAB script by inserting the `up.paths.root_folder` and `up.paths.save_folder` into the `setup_up` function.
   - _Note:_ The duration of recordings provided by the script can be adjusted with the `up.settings.req_durn` variable within the `setup_up` function.
3. Run the MATLAB script to download the required files from PhysioNet, import the data into MATLAB, and collate the data into a single MATLAB data file, ready for analysis.
 
