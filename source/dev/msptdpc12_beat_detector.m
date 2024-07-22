function [peaks, onsets] = msptdpc12_beat_detector(sig,fs)

% version: reduced lms scales (40 bpm)

options.use_reduced_lms_scales = true;

options.plaus_hr_bpm = [40, 200];

[peaks, onsets] = msptdpcref_beat_detector(sig,fs, options);

end
