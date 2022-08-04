# `AMPD_BEAT_DETECTOR` - AMPD PPG beat detector.
AMPD_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
using the 'Automatic Multiscale-based Peak Detection' beat detector

##  Inputs
+   sig : a vector of PPG values
    
+   fs  : the sampling frequency of the PPG in Hz
    
+   lpf_freq : (optional) the frequency (in Hz) of low-pass filtering applied to the PPG
    
##  Outputs
+   peaks : indices of detected pulse peaks
    
+   onsets : indices of detected pulse onsets
    
##  Reference
F. Scholkmann et al., 'An efficient algorithm for automatic peak detection in noisy periodic and quasi-periodic signals,' Algorithms, vol. 5, no. 4, pp. 588-603, 2012. <https://doi.org/10.3390/a5040588>

##  Author
Peter H. Charlton, University of Cambridge, February 2022.

##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Version
1.0

##  License - MIT
       Copyright (c) 2022 Peter H. Charlton
       Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
       The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
       THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
