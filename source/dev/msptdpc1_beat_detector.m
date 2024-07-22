function [peaks, onsets] = msptdpc1_beat_detector(sig,fs)

% version: pks only

options.find_pks = 1;
options.find_trs = 0;

[peaks, onsets] = msptdpcref_beat_detector(sig,fs, options);

end
