# Contributing to PPG-beats

Thanks for your interest in contributing!

The following guidelines are intended to provide an efficient way to contribute to the project. They are just guidelines, and if you think of better or additional ways to contribute then feel free to propose changes to this document in a pull request.

## Contributing to the software

Contributions to the software are most welcome. Here are some guidelines for common types of contributions.

### New PPG beat detection algorithms

We very much value contributions of new PPG beat detection algorithms. If you are considering contributing a new beat detection algorithm, then please:

1. Check to see whether this algorithm is already available in the toolbox by consulting the list of existing beat detectors [here](https://ppg-beats.readthedocs.io/en/latest/toolbox/ppg_beat_detectors/).
2. Prepare a submission:
   - **Code:** The algorithm should be provided in Matlab format, in a `.m` file, laid out as follows:
   ``` function [peaks, onsets] = <beat detector abbreviation>_beat_detector(sig,fs)
      
      <documentation>
      
      <code>
      
      end ```
      
      Note the following:
         - the main function takes the PPG signal (`sig`, a numerical column vector) and its sampling frequency (`fs`, a numerical variable) as inputs.
         - the main function is named "<beat detector abbreviation>_beat_detector".
         - the main function outputs two variables: `peaks` and `onsets`, which are column vectors containing the indices of detected pulse peaks / onsets in the signal.
         - additional functions can be used in the same file, and should appear below this main function.

### Specific improvements

### New feature suggestions

## Reporting issues or problems with the software

e.g. bugs

## Seeking support

## Improving the documentation

---

**Acknowledgment:** _This document is inspired in part by the [guidelines for contributing to Atom](https://github.com/atom/atom/blob/master/CONTRIBUTING.md).