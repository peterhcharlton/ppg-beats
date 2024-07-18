function [peaks, onsets] = msptdpc8_beat_detector(sig,fs)

% version: ds to 10 hz

options.do_ds = true;
options.ds_freq = 10;

[peaks, onsets] = msptdpcref_beat_detector(sig,fs, options);

end
