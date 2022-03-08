# Getting started with PPG-beats

A few pointers on how to quickly start using the toolbox.

---

## Requirements

The toolbox is run in Matlab, and requires the following add-on:

- **Signal Processing Toolbox**: This toolbox is used by both the performance assessment code, by several PPG beat detectors, and by an ECG beat detector. Therefore, it is required. If necessary, its use could be avoided if not pre-processing the PPG signal (as it is used several times by `preprocess_ppg_signal`), and if not using those beat detectors which require it.

In addition, the following add-ons are also required to enable the full functionality of the toolbox:

- **Statistics and Machine Learning Toolbox**: This toolbox is used by the `assess_multiple_datasets` code, and is only required when using this function to compare results across multiple datasets. It is also used by a couple of PPG beat detectors, although since these only use the `quantile` function from this toolbox, it wouldn't be too tricky to avoid the need for it in these beat detectors.
- **Econometrics Toolbox**: This toolbox is only used by one PPG beat detector.
- **Wavelet Toolbox**: This toolbox is only used by one or two PPG beat detectors.

_NB: You can obtain details of which functions use which toolboxes by running a Dependency Report in Matlab._

## Installation

The toolbox can be installed by:

1. Opening Matlab
2. Setting the current directory as the one where you want to save the toolbox, _e.g._
```cd C:/directoryname/```
3. Entering the following commands at the Matlab command window:

```
[old_path]=which('assess_beat_detectors'); if(~isempty(old_path)) rmpath(old_path(1:end-8)); end
toolbox_url='https://github.com/peterhcharlton/ppg-beats/archive/refs/heads/main.zip';
[filestr,status] = urlwrite(toolbox_url,'main.zip');
unzip('main.zip');
cd ppg-beats-main
addpath(genpath(pwd))
savepath
```
_NB: These instructions are adapted from those provided for the WFDB Toolbox [here](https://archive.physionet.org/physiotools/matlab/wfdb-app-matlab/)._

## Manual Download

The toolbox can be downloaded as a [ZIP file from GitHub](https://github.com/peterhcharlton/ppg-beats/archive/refs/heads/main.zip). It is available at:

- [GitHub](https://github.com/peterhcharlton/ppg-beats/)
- [Mathworks File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/107644-ppg-beats)
- [Zenodo](https://doi.org/10.5281/zenodo.6037646)

## Detecting beats in the PPG

The toolbox contains several PPG beat detectors, which are detailed [here](/toolbox/ppg_beat_detectors/).

This [PPG Beat Detection Tutorial](/tutorials/ppg_beat_detection/) provides instructions and sample code to quickly start using a PPG beat detection algorithm.

## Assessing the performance of PPG beat detectors

The toolbox contains tools to assess the performance of PPG beat detectors, which are detailed [here](/toolbox/performance_assessment/).

This [Performance Assessment Tutorial](/tutorials/performance_assessment/) provides instructions and sample code to quickly start using a PPG beat detection algorithm.
