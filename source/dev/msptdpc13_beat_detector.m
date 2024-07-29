function [peaks, onsets] = msptdpc13_beat_detector(sig,fs)

% version: optimal selection (pre 2024-07-01)

options.find_trs = true;
options.find_pks = false;
options.do_ds = true;
options.ds_freq = 20;
options.use_reduced_lms_scales = true;  % default 30 bpm

[peaks, onsets] = msptdpcref_beat_detector(sig,fs, options);

end
