# MIMIC 'PERform' Testing Dataset

---

This dataset contains ECG and PPG recordings of 10-minute duration, some of which were acquired from adults, and some of which were acquired neonates.

## Overview

 Item | Details 
 :--- | :--- 
 **Links** | [Dataset](https://doi.org/10.5281/zenodo.6807402), [Publication](http://peterhcharlton.github.io/publication/assess_ppg_beat_detectors/)
 **Signals** | PPG, ECG, resp
 **No. Subjs** | 200 (100 adults, 100 neonates)
 **Protocol** | Critically-ill patients during routine clinical care. Data were measured using a bedside monitor at 125 Hz. Data were extracted from the [MIMIC-III Waveform Database](https://physionet.org/content/mimic3wdb/1.0/).

## Downloading the dataset

Links to download the dataset in different formats are provided in this table:
 
Format | Links 
:--- | :--- 
Matlab | [all subjects](https://zenodo.org/record/6950488/files/mimic_perform_test_all_data.mat?download=1), [adult subjects](https://zenodo.org/record/6950488/files/mimic_perform_test_a_data.mat?download=1), [neonate subjects](https://zenodo.org/record/6950488/files/mimic_perform_test_n_data.mat?download=1)
CSV | [all subjects](https://zenodo.org/record/6950488/files/mimic_perform_test_all_csv.zip?download=1), [adult subjects](https://zenodo.org/record/6950488/files/mimic_perform_test_a_csv.zip?download=1), [neonate subjects](https://zenodo.org/record/6950488/files/mimic_perform_test_n_csv.zip?download=1)
WFDB | [all subjects](https://zenodo.org/record/6950488/files/mimic_perform_test_all_wfdb.zip?download=1), [adult subjects](https://zenodo.org/record/6950488/files/mimic_perform_test_a_wfdb.zip?download=1), [neonate subjects](https://zenodo.org/record/6950488/files/mimic_perform_test_n_wfdb.zip?download=1)

## Recreating the dataset

The dataset can be recreated from scratch by taking the following steps:

1. Use the following MATLAB script to import the data: [`collate_mimic_perform_train_test_datasets`](/functions/collate_mimic_perform_train_test_datasets/).
2. Modify the MATLAB script by inserting the `up.paths.root_folder` and `up.paths.save_folder` into the `setup_up` function.
   - _Note:_ The duration of recordings provided by the script can be adjusted with the `up.settings.req_durn` variable within the `setup_up` function.
   - _Note:_ The number of subjects in each group can be adjusted with the `up.settings.no_subjs_per_group` variable within the `setup_up` function.
3. Run the MATLAB script to download the required files from PhysioNet, import the data into MATLAB, and collate the data into a single MATLAB data file, ready for analysis.
4. (optionally) convert the dataset into WFDB and/or CSV formats using this MATLAB script: [`mimic_perform_dataset_converter`](/functions/mimic_perform_dataset_converter/)

