# `MSPTD_BEAT_DETECTOR` - MSPTD PPG beat detector.
MSPTD_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
using the 'Multi-Scale Peak and Trough Detection' beat detector

##  Inputs
+   sig : a vector of PPG values
    
+   fs  : the sampling frequency of the PPG in Hz
    
##  Outputs
+   peaks : indices of detected pulse peaks
    
+   troughs : indices of detected pulse troughs
    
##  Reference
S. M. Bishop and A. Ercole, 'Multi-scale peak and trough detection optimised for periodic and quasi-periodic neuroscience data,' in Intracranial Pressure and Neuromonitoring XVI. Acta Neurochirurgica Supplement, T. Heldt, Ed. Springer, 2018, vol. 126, pp. 189-195. <https://doi.org/10.1007/978-3-319-65798-1_39>

##  Author
Dr Steven Bishop - wrote the code ('new_peak_trough')
Peter H. Charlton - did very little, just wrote this wrapper and made slight edits (see 'new_peak_trough')

##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Version
0.1, and is still in development.

##  Licence
       ...TBC...
