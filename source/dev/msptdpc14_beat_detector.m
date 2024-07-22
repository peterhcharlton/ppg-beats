function [peaks, onsets] = msptdpc14_beat_detector(sig,fs)

% version: win len of 12 secs

options.win_len = 12;

[peaks, onsets] = msptdpcref_beat_detector(sig,fs, options);

end
