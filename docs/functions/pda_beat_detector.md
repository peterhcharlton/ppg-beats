# `PDA_BEAT_DETECTOR` - PDA PPG beat detector.
PDA_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
using the 'Peak Detection Algorithm' beat detector

##  Inputs
+   sig : a vector of PPG values
    
+   fs  : the sampling frequency of the PPG in Hz
    
##  Outputs
+   peaks : indices of detected pulse peaks
    
+   onsets : indices of detected pulse onsets
    
##  Reference
E. J. Arguello Prada and R. D. Serna Maldonado, 'A novel and low-complexity peak detection algorithm for heart rate estimation from low-amplitude photoplethysmographic (PPG) signals,' Journal of Medical Engineering and Technology, vol. 42, no. 8, pp. 569-577, 2018. <https://doi.org/10.1080/03091902.2019.1572237>

##  Author
+   Elisa Mejía Mejía - wrote the code: 'upslopes'
    
+   Peter H. Charlton - did very little, just wrote this wrapper
    
##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Version
1.0

##  License - MIT
       Copyright (c) 2022 Elisa Mejía Mejía and Peter H. Charlton
       Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
       The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
       THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
