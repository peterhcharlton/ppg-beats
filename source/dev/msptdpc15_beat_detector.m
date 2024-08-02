function [peaks, onsets] = msptdpc15_beat_detector(sig,fs)

% version: ds to 30 hz

options.do_ds = true;
options.ds_freq = 40;

[peaks, onsets] = msptdpcref_beat_detector(sig,fs, options);

end
