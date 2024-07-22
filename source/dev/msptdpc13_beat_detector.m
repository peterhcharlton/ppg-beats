function [peaks, onsets] = msptdpc13_beat_detector(sig,fs)

% version: optimal selection (pre 2024-07-01)

options.do_trs = true;
options.do_pks = false;
options.do_ds = true;
options.ds_freq = 20;
options.use_reduced_lms_scales = true;  % default 30 bpm

[peaks, onsets] = msptdpcref_beat_detector(sig,fs, options);

end
