# Performance Assessment

Resources to assess the performance of PPG beat detectors.

---

This page provides details of resources provided to assess the performance of PPG beat detectors.

## Datasets

Several [datasets](/datasets/summary) are provided on which to assess the performance of PPG beat detectors. Each dataset contains:

- PPG signals
- ECG signals (to obtain reference heartbeats)

## Code

**[`assess_beat_detectors.m`](/functions/assess_beat_detectors/)**: The main function for assessing the performance of PPG beat detectors. This can be used to assess the performance of beat detectors on a single dataset (as demonstrated in the tutorial [here](/tutorials/performance_assessment/)).

**[`assess_multiple_datasets.m`](/functions/assess_multiple_datasets/)**: This is a wrapper for `assess_beat_detectors.m`. It allows PPG beat detectors to be assessed across multiple datasets by: calling the main assessment function (`assess_beat_detectors.m`) for each dataset in turn; and then creating figures to display the performance on each dataset. [This tutorial](/tutorials/performance_assessment/) provides an example of how to use this function to assess performance across multiple datasets.