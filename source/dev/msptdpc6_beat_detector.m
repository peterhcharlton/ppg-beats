function [peaks, onsets] = msptdpc6_beat_detector(sig,fs)

% version: ds to 30 hz

options.do_ds = true;
options.ds_freq = 30;

[peaks, onsets] = msptdpcref_beat_detector(sig,fs, options);

end
