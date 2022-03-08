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



## Manual Download

...

## Detecting beats in the PPG

...

## Assessing the performance of PPG beat detectors

...
