# Contributing to PPG-beats

Thanks for your interest in contributing!

The following guidelines are intended to provide an efficient way to contribute to the project. They are just guidelines, and if you think of better or additional ways to contribute then feel free to propose changes to this document in a pull request.

## Contributing to the software

Contributions to the software are most welcome. Here are some guidelines for common types of contributions.

### New PPG beat detection algorithms

We very much value contributions of new PPG beat detection algorithms. If you are considering contributing a new beat detection algorithm, then please:

1. Check to see whether this algorithm is already available in the toolbox by consulting the list of existing beat detectors [here](https://ppg-beats.readthedocs.io/en/latest/toolbox/ppg_beat_detectors/).
2. Prepare a submission:
   
   **Code**
      
      The algorithm should be provided in Matlab format, in a `.m` file, laid out as follows:
      
      ```
      function [peaks, onsets] = <beat detector abbreviation>_beat_detector(sig,fs)
         
      <documentation>
         
      <code>
         
      end
      ```
      
      Note the following:
      
      - the main function takes the PPG signal (`sig`, a numerical column vector) and its sampling frequency (`fs`, a numerical variable) as inputs.
      - the main function is named "<beat detector abbreviation>_beat_detector".
      - the main function outputs two variables: `peaks` and `onsets`, which are column vectors containing the indices of detected pulse peaks / onsets in the signal.
      - additional functions can be used in the same file, and should appear below this main function.
      - the documentation in the file should include: (i) a reference to a publication describing the algorithm (if such a publication exists); (ii) the author's name; (iii) details of the license applied to the code, which should be an OSI-approved license (see [here](https://opensource.org/licenses/alphabetical)). The documentation will be automatically parsed to generate a webpage in the docs (e.g. [this webpage](https://ppg-beats.readthedocs.io/en/latest/functions/ampd_beat_detector/)), so please try to follow the documentation used in a template beat detector file ...
      - [This file](https://github.com/peterhcharlton/ppg-beats/blob/main/source/ampd_beat_detector.m) could be used as a template.
   
   **Description**
   
      Please prepare provide a short description of the algorithm (see [here](https://ppg-beats.readthedocs.io/en/latest/toolbox/ppg_beat_detectors/#automatic-multiscale-based-peak-detection) for an example) to be uploaded to the list of PPG beat detectors in the project documentation ([here](https://ppg-beats.readthedocs.io/en/latest/toolbox/ppg_beat_detectors/)).

3. Submit the new beat detector algorithm:
   
   **Code**
   
   Submit via a pull request.
   
   **Description**
   
   Open a new issue, referring to the submitted pull request, and providing the algorithm description.

### Specific improvements
   
You may well identify specific ways in which the code could be improved. Please submit revised code via a pull request.
   
_Improvements to beat detection algorithms:_ Note that most of the beat detection algorithms in this repository are either implementations of beat detectors which have been described in the literature, or implementations that have been previously provided by algorithm authors. Improvements which are intended to correct errors in these implementations will be considered. Improvements which change the design of the algorithm should instead be submitted following the approach described above for new beat detection algorithms.

### New feature suggestions
   
Please contribute new feature suggestions by opening a new issue.

## Reporting issues or problems with the software

Please report issues or problems with the software (e.g. bugs) by opening a new issue.

## Seeking support
   
If you would like support in using the code then please make your request by opening a new issue, and providing details of any problems you have encountered (ideally with an example which reproduces the problem). These will be addressed time-permitting.

## Improving the documentation
   
The documentation is written in markdown and stored in the [docs](https://github.com/peterhcharlton/ppg-beats/tree/main/docs) folder. Please submit improvements via a pull request.
   
_functions documentation_: Please note that all documentation in the [functions folder](https://github.com/peterhcharlton/ppg-beats/tree/main/docs/functions) is automatically generated from the comments within the matlab files in the [source folder](https://github.com/peterhcharlton/ppg-beats/tree/main/source). Therefore, any improvements to the functions documentation should be made by adjusting the matlab source files, rather than the automatically generated functions documentation.

---

**Acknowledgment:** _This document is inspired in part by the [guidelines for contributing to Atom](https://github.com/atom/atom/blob/master/CONTRIBUTING.md)._
