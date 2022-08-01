# PPG Beat Detectors

Algorithms to detect beats in photoplethysmogram (PPG) signals.

---

**Overview:** PPG-beats contains several algorithms to detect beats in the photoplethysmogram (PPG). This page provides an overview of these beat detectors. Follow the links for further details on each one.

**Accompanying tutorial:** See [this tutorial](./tutorials/ppg_beat_detection) for an example of how to use the beat detectors.

---

## Automatic Beat Detection

**Original publication:** Aboy M et al., An automatic beat detection algorithm for pressure signals. _IEEE Trans Biomed Eng_ 2005; 52: 1662-70. DOI: [10.1109/TBME.2005.855725](https://doi.org/10.1109/TBME.2005.855725)

**Description:** The PPG is strongly filtered to retain frequencies around an initial heart rate estimate, differentiated, and peaks are detected above the 75th percentile. Beats are identified as peaks in a weakly filtered PPG immediately following each peak identified in the differentiated signal.

**Link:** [abd_beat_detector](../../functions/abd_beat_detector)

**Licence:** GNU GPL Licence


## Automatic Multiscale-based Peak Detection

**Original publication:** F. Scholkmann et al., 'An efficient algorithm for
automatic peak detection in noisy periodic and quasi-periodic signals,'
_Algorithms_, vol. 5, no. 4, pp. 588-603, 2012. DOI: [10.3390/a5040588](https://doi.org/10.3390/a5040588)

**Description:** The PPG is detrended and segmented into 6s windows. A local maxima scalogram (LMS) is calculated: a matrix of random numbers, where the rows correspond to different scales (ranging from one sample to half the window duration), and the columns indicate PPG samples. The LMS values are set to zero when a PPG sample is higher than its neighbours at that particular scale. The LMS is truncated to only include scales smaller than the scale at which the most local maxima were identified. Beats are identified as samples which are deemed to be local maxima at all remaining scales.

**Link:** [ampd_beat_detector](../../functions/ampd_beat_detector)

**Licence:** MIT Licence


## Adaptive Threshold Beat Detector

**Original publication:** Shin HS et al., Adaptive threshold method for the peak detection of photoplethysmographic waveform. _Comput Biol Med_ 2009; 39: 1145-52. DOI: [10.1016/j.compbiomed.2009.10.006](https://doi.org/10.1016/j.compbiomed.2009.10.006)

**Description:** The PPG is bandpass filtered between 0.5 and 20 Hz. Troughs are identified as local minima which are below an adaptive threshold. The adaptive threshold increases from the value of the previous trough, at a rate related to the PPG amplitude. Any troughs occuring within a period of 0.6 times the previous inter-beat-interval are excluded. The 'Vmin' implementation of this beat detector was used, as it performed slightly better than the 'Vmax' implementation in initial testing.

**Link:** [atmax_beat_detector](../../functions/atmax_beat_detector) (see also: [atmin_beat_detector](../../functions/atmin_beat_detector))

**Licence:** MIT Licence


## Percentile Peak Detector

**Original publication:** C. Orphanidou et al., 'Signal-quality indices for the electrocardiogram and photoplethysmogram: derivation and applications to wireless monitoring,' _IEEE Journal of Biomedical and Health Informatics_, vol. 19, no. 3, pp. 832-8, 2015. DOI: [10.1109/JBHI.2014.2338351](https://doi.org/10.1109/JBHI.2014.2338351)

**Description:** In each 10 s PPG segment, beats are identified as peaks which are sufficiently close to (or above) the 90th percentile of the PPG signal, using adaptive filtering.

**Link:** [coppg_beat_detector](../../functions/coppg_beat_detector)

**Licence:** GNU GPL Licence


## Event-Related Moving Averages

**Original publication:** M. Elgendi et al., 'Systolic peak detection in acceleration
photoplethysmograms measured from emergency responders in tropical
conditions,' _PLoS ONE_, vol. 8, no. 10, pp. 1-11, 2013. DOI: [10.1371/journal.pone.0076585](https://doi.org/10.1371/journal.pone.0076585)

**Description:** The PPG is bandpass filtered between 0.5 and 8 Hz, rectified to eliminate values below zero, and squared. Two moving averages are calculated: (i) MApeak, a moving average of period 111 ms, emphasising systolic peaks; and (ii) MAbeat, a moving average of period 667 ms, emphasising individual beats. Beats are identified as maxima within periods lasting &ge; 111 ms where MApeak > MAbeat + &alpha; (where &alpha; is a threshold).

**Link:** [erma_beat_detector](../../functions/erma_beat_detector)

**Licence:** MIT Licence


## HeartPy

**Original publication:** P. van Gent et al., 'HeartPy: A novel heart rate algorithm for the
analysis of noisy signals,' _Transportation Research Part F: Traffic
Psychology and Behaviour_, vol. 66, pp. 368-378, 2019. DOI: [10.1016/j.trf.2019.09.015](https://doi.org/10.1016/j.trf.2019.09.015)

**Description:** The PPG is squared and normalised. Peaks are detected as maxima above a moving average (of period 0.75s). This is repeated for moving averages of different amplitudes, producing a set of peaks for each amplitude. The set of peaks which produces a plausible HR and the lowest variability in inter-beat intervals (IBIs) is selected as the set of beats. Beats which result in outlying IBIs are eliminated.

**Link:** [heartpy_beat_detector](../../functions/heartpy_beat_detector)

**Licence:** MIT Licence


## Incremental Merge Segmentation

**Original publication:** W. Karlen et al., 'Adaptive pulse segmentation and artifact detection in photoplethysmography for mobile applications,' in _Proc. IEEE EMBS_. IEEE, 2012, pp. 3131-4. DOI: [10.1109/EMBC.2012.6346628](https://doi.org/10.1109/EMBC.2012.6346628)

**Description:** Beats are detected at the end of continuous positive gradient segments (systolic upslopes) with an acceptable amplitude and duration, where the amplitude thresholds are adaptively calculated.

**Link:** [ims_beat_detector](../../functions/ims_beat_detector)

**Licence:** GNU GPL Licence


## Multi-Scale Peak and Trough Detection

**Original publication:** S. M. Bishop and A. Ercole, 'Multi-scale peak and trough detection optimised for periodic and quasi-periodic neuroscience data,' in _Intracranial Pressure and Neuromonitoring XVI. Acta Neurochirurgica Supplement_, T. Heldt, Ed. Springer, 2018, vol. 126, pp. 189-195. DOI: [10.1007/978-3-319-65798-1_39](https://doi.org/10.1007/978-3-319-65798-1_39)

**Description:** A modification of AMPD in which LMS matrices are calculated for both local maxima and local minima, so the algorithm detects both peaks and onsets. MSPTD also contains some optimisations to improve computational efficiency.

**Link:** [msptd_beat_detector](../../functions/msptd_beat_detector)

**Licence:** MIT Licence


## Peak Detection Algorithm

**Original publication:** E. J. Arguello Prada and R. D. Serna Maldonado, 'A novel and
low-complexity peak detection algorithm for heart rate estimation from low-amplitude photoplethysmographic (PPG) signals,' _Journal of Medical Engineering and Technology_, vol. 42, no. 8, pp. 569-577, 2018. DOI: [10.1080/03091902.2019.1572237](https://doi.org/10.1080/03091902.2019.1572237)

**Description:** Systolic peaks are identified as peaks which follow an upslope (i.e. period of positive gradient) lasting &ge; 60% of the duration of the upslope leading to the previously detected systolic peak.

**Link:** [pda_beat_detector](../../functions/pda_beat_detector)

**Licence:** MIT Licence


## Pulse Wave Delineator

**Original publication:** B. N. Li et al., 'On an automatic delineator for arterial blood pressure waveforms,' _Biomedical Signal Processing and Control_, vol. 5, no. 1, pp. 76-81, 2010. DOI: [10.1016/j.bspc.2009.06.002](https://doi.org/10.1016/j.bspc.2009.06.002)

**Description:** Pulse onsets and pulse peaks are identified from zero-crossing points in the first derivative of the PPG: onsets are identified as zero-crossing points before a maximal deflection, and peaks are identified as zero-crossing points immediately following maximal deflections.

**Link:** [pwd_beat_detector](../../functions/pwd_beat_detector)

**Licence:** The 2-Clause BSD License


## PPG Pulses Detector

**Original publication:** J. Lazaro et al., 'Pulse rate variability analysis for discrimination of sleep-apnea-related decreases in the amplitude fluctuations of pulse
photoplethysmographic signal in children,' _IEEE Journal of Biomedical and Health Informatics_, vol. 18, no. 1, pp. 240-246, 2014. DOI: [10.1109/JBHI.2013.2267096](https://doi.org/10.1109/JBHI.2013.2267096)

**Description:** Peaks are identified in the differentiated PPG using an adaptive filter set to the amplitude of the previous peak, and decreases for a period after that peak at a rate dependent on previous inter-beat intervals. Beats are identified as maxima in the PPG within 300ms of each peak in the differentiated PPG.

**Link:** [ppgpulses_beat_detector](../../functions/ppgpulses_beat_detector)

**Licence:** GNU GPL Licence


## Adapted Onset Detector

**Original publication:** A. N. Vest et al., 'An open source benchmarked toolbox for cardiovascular waveform and interval analysis,' _Physiological Measurement_, vol. 39, no. 10, 2018. DOI: [10.1088/1361-6579/aae021](https://doi.org/10.1088/1361-6579/aae021)

**Description:** Systolic upslopes are detected from a signal generated with a slope sum function, which sums the magnitudes of the PPG upslopes in the previous 0.17 s. Adaptive thresholding is used to identify systolic upslopes in this signal. The 'qppgfast' implementation of this beat detector was used, after testing showed it performed similarly to the original 'qppg' implementation.

**Links:** [qppg_beat_detector](../../functions/qppg_beat_detector) (original version); [qppgfast_beat_detector](../../functions/qppgfast_beat_detector) (faster version);

**Licences:** GNU GPL Licence (and also BSD 3-Clause License for 'qppgfast_beat_detector')


## Symmetric Projection Attractor Reconstruction Detector

**Original publication:**  n/a

**Description:** The PPG is segmented into 20 s windows and time delay coordinates are used to represent it in 7-dimensional phase space with the time delay set to one seventh of the average inter-beat interval. The Symmetric Projection Attractor Reconstruction method is then used to construct an appropriate 2-dimensional projection of the phase space. Beats are identified as times at which the orbit crosses the x-axis. This implementation uses information from previous windows to inform beat detections in the current window.

**Link:** [spar_beat_detector](../../functions/spar_beat_detector)

**Licence:** _TBC_


## Stationary Wavelet Transform Beat Detector

**Original publication:**  S. Vadrevu and M. Sabarimalai Manikandan, 'A robust pulse onset and peak detection method for automated PPG signal analysis system,' _IEEE Trans Instrum Meas_, vol.68, no.3, pp.807-817, 2019. DOI: [10.1109/TIM.2018.2857878](https://doi.org/10.1109/TIM.2018.2857878)

**Description:** The PPG is decomposed using the Stationary Wavelet Transform. Multi-scale sum and products of selected detail subbands are calculated to emphasise systolic upslopes. An envelope is then extracted by: adaptive thresholding to reduce the influence of noise; calculating the Shannon entropy; and smoothing the result. Finally, beats are identified in the envelope using a Gaussian derivative filter.

**Link:** [swt_beat_detector](../../functions/swt_beat_detector)

**Licence:** MIT Licence


## Wavelet Foot Delineation

**Original publication:** N. J. Conn and D. A. Borkholder, 'Wavelet based photoplethysmogram foot delineation for heart rate variability applications,' in _IEEE Signal Processing in Medicine and Biology Symposium_. IEEE, 2013. DOI: [10.1109/SPMB.2013.6736782](https://doi.org/10.1109/SPMB.2013.6736782)

**Description:** The PPG is bandpass filtered between 0.5 and 8 Hz, and interpolated to 250 Hz. It is decomposed using a wavelet transform, retaining the fifth wavelet scale for analysis. This signal is rectified and squared to eliminate values below zero. Regions containing beats are identified as those where the signal exceeds a low-pass filtered version of the signal. The timing of the beat within each region is identified as the first zero-crossing of the third derivative, or failing that, the maximum in the second derivative.

**Link:** [wfd_beat_detector](../../functions/wfd_beat_detector)

**Licence:** MIT Licence

---

_**Source:** details of the peak detectors are taken from [DOI: 10.1088/1361-6579/ac826d](http://peterhcharlton.github.io/publication/assess_ppg_beat_detectors/) under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)._

---