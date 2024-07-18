# Designing a Beat Detector

The tutorials on this page demonstrate the processes undertaken to design MSPTDfast, an efficient PPG beat detection algorithm.

## Designing MSPTDfast (v.1)

This tutorial demonstrates the processes undertaken to design MSPTDfast (v.1), which was designed using a single dataset.

- Install the PPG-beats toolbox by following the instructions [here](https://ppg-beats.readthedocs.io/en/latest/toolbox/getting_started/).
- Download the PPG-DaLiA dataset in Matlab format from [here](https://drive.google.com/uc?export=download&id=1FWjLR1-BaoQQtezityfvPQs53_d-xZTX).
- Use the `assess_multiple_datasets.m` script to run the analysis.

- During this process a new folder will have been created called `proc_data_ppg_dalia_lunch_break`. Within this folder you will find files storing the analysis steps, including the file containing the results: the `ppg_detect_stats.mat` file.
- Make a copy of this `ppg_detect_stats.mat` file and rename it as `ppg_detect_stats_lunchbreak.mat`. Note down the location of this file.
- Finally, analyse the performance of the different algorithm configuration options by running the `msptdfast_cinc_analysis.m` script.

This tutorial is demonstrated in the following video:

<iframe width="560" height="315" src="https://www.youtube.com/embed/jd5fT_HvRLo" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

