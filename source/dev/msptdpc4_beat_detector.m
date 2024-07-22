function [peaks, onsets] = msptdpc4_beat_detector(sig,fs)

% version: reduced lms scales (30 bpm)

options.use_reduced_lms_scales = true;

[peaks, onsets] = msptdpcref_beat_detector(sig,fs, options);

end
