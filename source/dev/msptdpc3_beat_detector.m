function [peaks, onsets] = msptdpc3_beat_detector(sig,fs)

% version: vectorisation 

options.lms_calc_method = 2;

[peaks, onsets] = msptdpcref_beat_detector(sig,fs, options);

end
