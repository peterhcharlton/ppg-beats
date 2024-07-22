# Designing a Beat Detector

The tutorials on this page demonstrate the processes undertaken to design MSPTDfast, an efficient PPG beat detection algorithm.

## Designing MSPTDfast (v.1)

This tutorial demonstrates the processes undertaken to design MSPTDfast (v.1), which was designed using a single dataset.

- Install the PPG-beats toolbox by following the instructions [here](https://ppg-beats.readthedocs.io/en/latest/toolbox/getting_started/). _NB: the [v.2.0](https://github.com/peterhcharlton/ppg-beats/releases/tag/v.2.0) release of the toolbox was used for the reported analysis, so should be used to replicate the analysis exactly._
- Download the PPG-DaLiA dataset in Matlab format from [here](https://zenodo.org/records/12793711/files/ppg_dalia_lunch_break_data.mat?download=1).
- Use the `assess_multiple_datasets.m` script to run the analysis.
- During this process a new folder will have been created called `proc_data_ppg_dalia_lunch_break`. Within this folder you will find files storing the analysis steps, including the file containing the results: the `ppg_detect_stats.mat` file. Note down the location of this file.
- Finally, analyse the performance of the different algorithm configuration options by running the `[msptdfast_cinc_analysis.m](https://raw.githubusercontent.com/peterhcharlton/ppg-beats/main/source/publication_specific_scripts/msptdfast_cinc_analysis_post_20240701.m)` script.

This tutorial is demonstrated in the following video:

<iframe width="560" height="315" src="https://www.youtube.com/embed/jd5fT_HvRLo" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

