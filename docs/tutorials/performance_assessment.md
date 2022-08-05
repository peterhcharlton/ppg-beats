# Performance Assessment

The tutorials on this page demonstrate how to assess the performance of PPG beat detectors using publicly available datasets.

## Assessing performance on a single dataset

This tutorial demonstrates how to assess the performance of beat detectors on a single dataset.

- Install the PPG-beats toolbox by following the instructions [here](https://ppg-beats.readthedocs.io/en/latest/toolbox/getting_started/).
- Download the truncated version of the MIMIC PERform Training Dataset from [here](https://zenodo.org/record/6967256/files/MIMIC_PERform_truncated_train_all_data.mat?download=1).
- Use the `assess_beat_detectors.m` script to run the analysis using the following commands:
```matlab
options.dataset_file = 'filepath'; % NB: Here, instead of 'filepath' insert the path of the `MIMIC_PERform_truncated_train_all_data.mat` file on your computer, _e.g._ `/Users/petercharlton/Downloads/MIMIC_PERform_truncated_train_all_data.mat`
assess_beat_detectors('MIMIC_PERform_truncated_train_all_data', options)
```

- During this process a new folder will have been created called `proc_data_MIMIC_PERform_truncated_train_all_data`. Within this folder you will find files storing the analysis steps, including the file containing the results: the `ppg_detect_stats.mat` file. Open this file in Matlab, and you will see it contains a variable called `res` (for results).
- Have a look at the table containing the results by typing `res.noQual.txt` at the command window. This provides the median and inter-quartile range for a set of statistics describing the performance of each beat detector (each 'strategy' in this table).

## Assessing performance on multiple datasets

_Tutorial to be inserted, using `assess_multiple_datasets.m`_