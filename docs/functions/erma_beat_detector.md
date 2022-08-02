# `ERMA_BEAT_DETECTOR` - ERMA PPG beat detector.
ERMA_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
using the 'Event-Related Moving Averages' beat detector

##  Inputs
+   sig : a vector of PPG values
    
+   fs  : the sampling frequency of the PPG in Hz
    
##  Outputs
+   peaks : indices of detected pulse peaks
    
+   onsets : indices of detected pulse onsets
    
##  Reference
M. Elgendi et al., 'Systolic peak detection in acceleration photoplethysmograms measured from emergency responders in tropical conditions,' PLoS ONE, vol. 8, no. 10, pp. 1-11, 2013. <https://doi.org/10.1371/journal.pone.0076585>

##  Author
+   Elisa Mejía Mejía - wrote the code: 'erma'
    
+   Peter H. Charlton - did very little, just wrote this wrapper and made slight edits (see 'erma')
    
##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Version
1.0

##  MIT License
       Copyright (c) 2022 Elisa Mejía Mejía and Peter H. Charlton
       Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
       The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
       THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
