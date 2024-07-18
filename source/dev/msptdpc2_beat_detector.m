function [peaks, onsets] = msptdpc2_beat_detector(sig,fs)

% version: trs only

options.find_pks = 0;
options.find_trs = 1;

[peaks, onsets] = msptdpcref_beat_detector(sig,fs, options);

end
