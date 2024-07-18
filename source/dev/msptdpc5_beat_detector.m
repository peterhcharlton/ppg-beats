function [peaks, onsets] = msptdpc5_beat_detector(sig,fs)

% version: default

options = struct;
[peaks, onsets] = msptdpcref_beat_detector(sig,fs, options);

end
