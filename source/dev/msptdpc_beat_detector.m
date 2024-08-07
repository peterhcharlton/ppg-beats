function [peaks, onsets] = msptdpc_beat_detector(sig,fs)
% MSPTD_BEAT_DETECTOR  MSPTD PPG beat detector.
%   MSPTD_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
%   using the 'Multi-Scale Peak and Trough Detection' beat detector
%   
%   # Inputs
%   
%   * sig : a vector of PPG values
%   * fs  : the sampling frequency of the PPG in Hz
%   
%   # Outputs
%   * peaks : indices of detected pulse peaks
%   * onsets : indices of detected pulse troughs (i.e. onsets)
%   
%   # Reference
%   S. M. Bishop and A. Ercole, 'Multi-scale peak and trough detection optimised for periodic and quasi-periodic neuroscience data,' in Intracranial Pressure and Neuromonitoring XVI. Acta Neurochirurgica Supplement, T. Heldt, Ed. Springer, 2018, vol. 126, pp. 189-195. <https://doi.org/10.1007/978-3-319-65798-1_39>
%   
%   # Author
%   Peter H. Charlton
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

%% detect peaks and onsets in each window

peaks = []; onsets = [];
for win_no = 1 : length(win_starts)

    % - extract this window's data
    win_sig = sig(win_starts(win_no):win_ends(win_no));
    
    % - downsample signal
    if do_ds
        rel_sig = downsample(win_sig, ds_factor);
        rel_fs = ds_fs;
    else
        rel_sig = win_sig;
        rel_fs = fs;
    end
    
    % - detect peaks and onsets
    [p,t] = detect_peaks_and_onsets_using_msptd(win_sig,rel_fs);
    
    % - resample peaks
    if do_ds
        p = p*ds_factor;
        t = t*ds_factor;
    end
    
    % - correct peak indices by finding highest point within tolerance either side of detected peaks
    tol = ceil(fs*0.05);
    for pk_no = 1 : length(p)
        [~, temp] = max(rel_sig( (p(pk_no) - tol) : (p(pk_no) + tol) ));
        p(pk_no) = p(pk_no) - tol + temp - 1;
    end
    
    % - correct onset indices by finding highest point within tolerance either side of detected onsets
    for onset_no = 1 : length(t)
        [~, temp] = min(rel_sig( (t(onset_no) - tol) : (t(onset_no) + tol) ));
        t(onset_no) = t(onset_no) - tol + temp - 1;
    end
    
    % - store peaks and onsets
    win_peaks = p + win_starts(win_no) -1;
    peaks = [peaks; win_peaks];
    win_onsets = t + win_starts(win_no) -1;
    onsets = [onsets; win_onsets];
    
    clear temp
    
end

% tidy up detected peaks and onsets (by ordering them and only retaining unique ones)
peaks = unique(peaks(:));
onsets = unique(onsets(:));

%% correct peak and onset indices

% (i) to account for downsampling
if do_ds
    % - peaks
    peaks = peaks*ds_factor;
    % find highest point within tolerance either side of detected peaks
    for pk_no = 1 : length(peaks)
        curr_peak = peaks(pk_no);
        tol_start = curr_peak - ds_factor;
        tol_end = curr_peak + ds_factor;
        [~, temp] = max(sig(tol_start:tol_end));
        peaks(pk_no) = curr_peak - ds_factor + temp - 1;
        clear temp curr_peak tol_start tol_end
    end
    % - onsets
    onsets = onsets*ds_factor;
    % find highest point within tolerance either side of detected onsets
    for onset_no = 1 : length(onsets)
        curr_onset = onsets(onset_no);
        tol_start = curr_onset - ds_factor;
        tol_end = curr_onset + ds_factor;
        [~, temp] = min(sig(tol_start:tol_end));
        onsets(onset_no) = curr_onset - ds_factor + temp - 1;
        clear temp curr_onset tol_start tol_end
    end
end

% (ii) tidy up
peaks = tidy_beats(peaks);
onsets = tidy_beats(onsets);

%% check
do_check = 0;
if do_check
    plot(sig), hold on, 
    plot(peaks, sig(peaks), '*r')
    plot(onsets, sig(onsets), '*k')
    close all
end

end

function [p,t] = detect_peaks_and_onsets_using_msptd(x,fs)

% Setup
N = length(x); % length of signal
L = ceil(N/2)-1; % max window length

% Step 0: (added by PC) don't calculate scales outside the range of plausible HRs
plaus_hr_hz = [30, 200]./60; % in Hz
init_scales = 1:L;
durn_signal = length(x)/fs;
init_scales_fs = (L./init_scales)/durn_signal;
%init_scales_inc_log = init_scales_fs >= plaus_hr_hz(1) & init_scales_fs <= plaus_hr_hz(2);
init_scales_inc_log = init_scales_fs >= plaus_hr_hz(1);
max_scale = find(init_scales_inc_log,1,'last');
%init_scales_inc_log = true(size(init_scales_inc_log));

% Step 1: calculate local maxima and local minima scalograms

% - detrend
x = detrend(x); % this removes the best-fit straight line

% - initialise LMS matrices
% m_max = false(L,N);
% m_min = false(L,N);
m_max = false(max_scale,N);
m_min = false(max_scale,N);

% - populate LMS matrices
for k = 1:max_scale % scalogram scales
    for i = (k+2):(N-k+1)
        if x(i-1) > x(i-k-1) && x(i-1) > x(i+k-1)
            m_max(k,i-1) = true;
        end
        if x(i-1) < x(i-k-1) && x(i-1) < x(i+k-1)
            m_min(k,i-1) = true;
        end
    end
end

% Step 2: find the scale with the most local maxima (or local minima)
% - row-wise summation (i.e. sum each row)
gamma_max = sum(m_max,2); % the ",2" option makes it row-wise
gamma_min = sum(m_min,2); % the ",2" option makes it row-wise
% - find scale with the most local maxima (or local minima)
[~, lambda_max] = max(gamma_max);
[~, lambda_min] = max(gamma_min);

% Step 3: Use lambda to remove all elements of m for which k>lambda
first_scale_to_include = find(init_scales_inc_log,1);
m_max = m_max(first_scale_to_include:lambda_max,:);
m_min = m_min(first_scale_to_include:lambda_min,:);

% Step 4: Find peaks
% - column-wise summation
m_max_sum = sum(~m_max);
m_min_sum = sum(~m_min);
p = find(~m_max_sum); p = p(:);
t = find(~m_min_sum); t = t(:);
end
