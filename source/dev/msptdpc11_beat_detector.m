function [peaks, onsets] = msptdpc11_beat_detector(sig,fs)

% version: win len of 10 secs

options.win_len = 10;

[peaks, onsets] = msptdpcref_beat_detector(sig,fs, options);

end
