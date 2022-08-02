# `PWD_BEAT_DETECTOR` - PWD PPG beat detector.
PWD_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
using the 'Pulse Wave Delineator' beat detector

##  Inputs
+   sig : a vector of PPG values
    
+   fs  : the sampling frequency of the PPG in Hz
    
##  Outputs
+   peaks : indices of detected pulse peaks
    
+   onsets : indices of detected pulse onsets
    
##  Reference
B. N. Li et al., 'On an automatic delineator for arterial blood pressure waveforms,' Biomedical Signal Processing and Control, vol. 5, no. 1, pp. 76-81, 2010. <https://doi.org/10.1016/j.bspc.2009.06.002>

##  Author
+   LI Bing Nan - created the original code
    
+   Peter H. Charlton - did very little, just wrote this wrapper
    
##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Version
1.0

##  Source
<https://www.mathworks.com/matlabcentral/fileexchange/29484-pulse-waveform-delineator>

##  Licence
       The 2-Clause BSD License (see the copyright notice, list of conditions, and disclaimer below)
