function [peaks, onsets] = msptdpc9_beat_detector(sig,fs)

% version: win len of 4 secs

options.win_len = 4;

[peaks, onsets] = msptdpcref_beat_detector(sig,fs, options);

end
