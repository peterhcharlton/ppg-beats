# `PPGPULSES_BEAT_DETECTOR` - PPGPULSES PPG beat detector.
PPGPULSES_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
using the 'PPG Pulses' beat detector

##  Inputs
+   sig : a vector of PPG values
    
+   fs  : the sampling frequency of the PPG in Hz
    
##  Outputs
+   peaks : indices of detected pulse peaks
    
+   onsets : indices of detected pulse onsets
    
##  Reference
J. Lazaro et al., 'Pulse rate variability analysis for discrimination of sleep-apnea-related decreases in the amplitude fluctuations of pulse photoplethysmographic signal in children,' IEEE Journal of Biomedical and Health Informatics, vol. 18, no. 1, pp. 240-246, 2014. <https://doi.org/10.1109/JBHI.2013.2267096>

##  Author
+   Jesus Lazaro - created the original code
    
+   Mariano Llamedo Soria - adapted it to the ecg-kit project
    
+   Peter H. Charlton - did very little, just wrote this wrapper
    
##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Version
1.0

##  Source
This function was copied from the 'ecg-kit' repository by 'marianux' at:
<https://github.com/marianux/ecg-kit> , where it is licensed under the 
GNU License v2.0. See also: <http://marianux.github.io/ecg-kit/>
This particular function was available at: <https://github.com/marianux/ecg-kit/blob/master/common/ppg/PPG_pulses_detector.m>

##  Licence
       This is available under the GNU General Public License v2.0: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
