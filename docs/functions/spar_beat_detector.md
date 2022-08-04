# `SPAR_BEAT_DETECTOR` - SPAR PPG beat detector.
SPAR_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
using the 'Symmetric Projection Attractor Reconstruction' beat detector

##  Inputs
+   sig : a vector of PPG values
    
+   fs  : the sampling frequency of the PPG in Hz
    
##  Outputs
+   peaks : indices of detected pulse peaks
    
+   onsets : indices of detected pulse onsets
    
##  Reference
C. Pettit and P. Aston, 'Photoplethysmogram (PPG) Beat Detection Using Symmetric Projection Attractor Reconstruction,' [in preparation].

##  Author
+   Callum Pettit and Philip Aston - wrote the code: 'SPARcycleTimesPPG' and 'find_acl'
    
+   Peter H. Charlton - modified the code slightly (as indicated), but mainly just wrote this wrapper
    
##  Documentation
<https://ppg-beats.readthedocs.io/>

##  Version
1.0

##  Licence
       SPAR PPG Beat Detector as [downloaded/received] is Copyright ©, 2021-2022 University of Surrey.
    
       All Rights Reserved - You are allowed to alter and use this code only for your own non-commercial research purposes and this excludes commercial use of any kind. You are not permitted to copy and distribute this code to anyone under any circumstances. All rights remain with the Copyright holders.
    
       If you wish to copy, distribute or make any commercial use of the Software please contact Surrey’s IP & Licensing team at techtransfer@surrey.ac.uk to request an appropriate licence.
    
       Please contact Philip Aston (P.Aston@surrey.ac.uk) with any technical questions or collaboration requests.
