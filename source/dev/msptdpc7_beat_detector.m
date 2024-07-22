function [peaks, onsets] = msptdpc7_beat_detector(sig,fs)

% version: ds to 20 hz

options.do_ds = true;
options.ds_freq = 20;

[peaks, onsets] = msptdpcref_beat_detector(sig,fs, options);

end
