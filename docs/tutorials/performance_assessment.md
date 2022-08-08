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

This tutorial is demonstrated in the following video:

<iframe width="560" height="315" src="https://www.youtube.com/embed/jd5fT_HvRLo" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>


## Assessing performance on multiple datasets

This tutorial demonstrates how to assess the performance of beat detectors across multiple datasets.

- Install the PPG-beats toolbox by following the instructions [here](https://ppg-beats.readthedocs.io/en/latest/toolbox/getting_started/).
- Download the required datasets:
    - [MIMIC PERform Truncated Training dataset](https://zenodo.org/record/6967256/files/MIMIC_PERform_truncated_train_all_data.mat?download=1)
    - [MIMIC PERform Truncated Testing dataset](https://zenodo.org/record/6973963/files/MIMIC_PERform_truncated_test_all_data.mat?download=1)
    - MIMIC PERform Truncated AF dataset (both [AF subjects](https://zenodo.org/record/6973963/files/MIMIC_PERform_truncated_af_data.mat?download=1) and [non-AF subjects](https://zenodo.org/record/6973963/files/MIMIC_PERform_truncated_non_af_data.mat?download=1)).
- Open the `assess_multiple_datasets.m` script (one of the scripts within the **PPG-beats** toolbox) in Matlab. Edit it as follows:
    - Specify the path of the Matlab file for each dataset in the `specify_path_of_dataset_file` function within the `assess_multiple_datasets.m` script.
    - Specify the datasets to be analysed. This can be quickly done by enabling the `do_tutorial` option within the `setup_universal_params` function within the `assess_multiple_datasets.m` script.
    - Specify the folder in which to save results figures by modifying the `up.paths.plots_root_folder` variable within the `setup_universal_params` function within the `assess_multiple_datasets.m` script.
- Run the `assess_multiple_datasets.m` script to assess the performance of PPG beat detectors across all of the datasets. This will:
    - Analyse each dataset in turn
    - Generate results (in the command window) and results figures (saved in the specified folder).

This tutorial is demonstrated in the following video:

<iframe width="560" height="315" src="https://www.youtube.com/embed/Li7aiqPSSfc" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>