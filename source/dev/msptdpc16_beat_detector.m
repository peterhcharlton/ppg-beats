function [peaks, onsets] = msptdpc16_beat_detector(sig,fs)

% version: ds to 30 hz

options.do_ds = true;
options.ds_freq = 50;

[peaks, onsets] = msptdpcref_beat_detector(sig,fs, options);

end
