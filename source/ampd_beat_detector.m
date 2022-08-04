function [peaks, onsets] = ampd_beat_detector(sig,fs)
% AMPD_BEAT_DETECTOR  AMPD PPG beat detector.
%   AMPD_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
%   using the 'Automatic Multiscale-based Peak Detection' beat detector
%   
%   # Inputs
%   
%   * sig : a vector of PPG values
%   * fs  : the sampling frequency of the PPG in Hz
%   * lpf_freq : (optional) the frequency (in Hz) of low-pass filtering applied to the PPG
%   
%   # Outputs
%   * peaks : indices of detected pulse peaks
%   * onsets : indices of detected pulse onsets
%   
%   # Reference
%   F. Scholkmann et al., 'An efficient algorithm for automatic peak detection in noisy periodic and quasi-periodic signals,' Algorithms, vol. 5, no. 4, pp. 588-603, 2012. <https://doi.org/10.3390/a5040588>
%   
%   # Author
%   Peter H. Charlton, University of Cambridge, February 2022.
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Version
%   1.0
%   
%   # License - MIT
%      Copyright (c) 2022 Peter H. Charlton
%      Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
%      The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
%      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%% Window signal
% split into overlapping 6 s windows
win_len = 6; % in secs
overlap = 0.2; % proportion of overlap between consecutive windows
no_samps_in_win = win_len*fs;
if length(sig) <= no_samps_in_win
    win_starts = 1;
    win_ends = length(sig);
else
    win_offset = round(no_samps_in_win*(1-overlap));
    win_starts = 1:win_offset:length(sig)-no_samps_in_win;
    win_ends = win_starts + no_samps_in_win;
    if win_ends(end) < length(sig)
        win_starts(end+1) = length(sig) - no_samps_in_win;
        win_ends(end+1) = length(sig);
    end % this ensures that the windows include the entire signal duration
end

%% Downsample signal

% Set up downsampling if the sampling frequency is particularly high
do_ds = 0;
if do_ds
    min_fs = 2*10;
    if fs > min_fs
        ds_factor = floor(fs/min_fs);
        ds_fs = fs/floor(fs/min_fs);
    else
        do_ds = 0;
    end
end

%% detect peaks in each window

peaks = [];
for win_no = 1 : length(win_starts)

    % - extract this window's data
    win_sig = sig(win_starts(win_no):win_ends(win_no));
    
    % - downsample signal
    if do_ds
        rel_sig = downsample(win_sig, ds_factor);
    else
        rel_sig = win_sig;
    end
    
    % - detect peaks
    p = detect_peaks_using_ampd(win_sig);
    
    % - resample peaks
    if do_ds
        p = p*ds_factor;
    end
    
    % - store peaks
    win_peaks = p + win_starts(win_no) -1;
    peaks = [peaks; win_peaks];
end

% tidy up detected peaks and onsets (by ordering them and only retaining unique ones)
peaks = unique(peaks(:));

%% correct peak indices

% (i) peaks always seem to be one index too high
peaks = peaks-1;

% (ii) to account for downsampling
if do_ds
    peaks = peaks*ds_factor;
    % find highest point within tolerance either side of detected peaks
    for pk_no = 1 : length(peaks)
        curr_peak = peaks(pk_no);
        tol_start = curr_peak - ds_factor;
        tol_end = curr_peak + ds_factor;
        [~, temp] = max(sig(tol_start:tol_end));
        peaks(pk_no) = curr_peak - ds_factor + temp;
        clear temp curr_peak tol_start tol_end

    end
end

% (iii) tidy up
peaks = tidy_beats(peaks);

%% calculate onset indices
onsets = pulse_onsets_from_peaks(sig, peaks);

%% check
do_check = 0;
if do_check
    plot(sig), hold on, 
    plot(peaks, sig(peaks), '*r')
    plot(onsets, sig(onsets), '*k')
    close all
end

end

function p = detect_peaks_using_ampd(x)
N = length(x); % length of signal
L = ceil(N/2)-1; % max window length
alpha = 1;  % constant factor

% Step 1: calculate local maxima scalogram (LMS)

% - detrend
x = detrend(x); % this removes the best-fit straight line

% - initialise LMS matrix
m = alpha + rand(L,N); % rand is uniform across interval (0,1)

% - populate LMS matrix
for k = 1:L % local maxima scalogram scales
    for i = (k+2):(N-k+1)
        if x(i-1) > x(i-k-1) && x(i-1) > x(i+k-1)
            m(k,i) = 0;
        end
    end
end

% Step 2: find the scale with the most local maxima
% - row-wise summation (i.e. sum each row)
gamma = sum(m,2); % the ",2" option makes it row-wise
% - find scale with the most local maxima (which is the lowest value of gamma as maxima are set to zero).
[~, lambda] = min(gamma);

% Step 3: Use lambda to remove all elements of m for which k>lambda
m = m(1:lambda,:);

% Step 4: Find peaks
sigma = std(m); % column-wise standard deviation
p = find(sigma==0);
p = p(:);

end
