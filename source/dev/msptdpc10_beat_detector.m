function [peaks, onsets] = msptdpc10_beat_detector(sig,fs)

% version: win len of 6 secs

options.win_len = 6;

[peaks, onsets] = msptdpcref_beat_detector(sig,fs, options);

end
