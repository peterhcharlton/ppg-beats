# `HEARTPY_BEAT_DETECTOR` - HEARTPY PPG beat detector.
HEARTPY_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
using the 'HeartPy' beat detector

##  Inputs
+   sig : a vector of PPG values
    
+   fs  : the sampling frequency of the PPG in Hz
    
##  Outputs
+   peaks : indices of detected pulse peaks
    
+   onsets : indices of detected pulse onsets
    
##  Reference
P. van Gent et al., 'HeartPy: A novel heart rate algorithm for the analysis of noisy signals,' Transportation Research Part F: Traffic Psychology and Behaviour, vol. 66, pp. 368-378, 2019. https://doi.org/10.1016/j.trf.2019.09.015>

##  Author
Peter H. Charlton, University of Cambridge, February 2022.

##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Version
1.0

##  Source
This is a Matlab implementation of the original Python algorithm, described in:
van Gent P. et al., Analysing noisy driver physiology real-time using off-the-shelf sensors: Heart rate analysis software from the taking the fast lane project. J. Open Res. Softw. 2019, 7, <https://doi.org/10.5334/jors.241>
It was created based on v.1.2.5 of HeartPy downloaded on 7-June-2021 from: <https://github.com/paulvangentcom/heartrate_analysis_python> 
NB: The original HeartPy algorithm contains several extra features (such as artifact detection), which are not included in this code.

##  MIT License
       Original:
       MIT License
       Copyright (c) 2021 Paul van Gent
       Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
       The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
       THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    
       The code has since been implemented in Matlab format, and modified, by Peter H. Charlton.
