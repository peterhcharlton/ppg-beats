# ECG Beat Detection Algorithms

---

The performance of PPG beat detection algorithms in **PPG-beats** can be assessed by comparing the beat detections against reference beats obtained from simultaneous electrocardiogram (ECG) signals. To facilitate this, **PPG-beats** contains algorithms to detect beats in the ECG, and assess the quality of these beat detections.

## _gqrs_ Beat Detector

**Original publication:** The methodology used by this beat detector isn't described in full in a publication. However, I've found [this blog post](https://c4science.ch/source/Event_based_gQRS/browse/master/python_gQRS/gQRS-lvlCrossingsampling.ipynb) helpful for understanding how it works. I would also recommend reading [this publication](https://doi.org/10.1088/0967-3334/36/8/1665) and [this assessment of ECG beat detectors](https://ieeexplore.ieee.org/document/7043053).

**Description:** 

**Link:** [c-code](https://www.physionet.org/physiotools/wag/gqrs-1.htm), [MATLAB wrapper](https://archive.physionet.org/physiotools/matlab/wfdb-app-matlab/html/gqrs.m)

## _jqrs_ Beat Detector

**Original publications:** 

1. J. Behar et al., 'A comparison of single channel fetal ecg extraction methods,' _Ann Biomed Eng_ 2014; 42: 1340-53. DOI: [10.1007/s10439-014-0993-9](https://doi.org/10.1007/s10439-014-0993-9)
2. A.E.W. Johnson _et al._, 'R-peak estimation using multimodal lead switching,' _Proc CinC._, 2014. p. 281-4. [Link](https://ieeexplore.ieee.org/document/7043034)

**Description:** 

**Link:** 

## Quality Assessment

