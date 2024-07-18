function [peaks, onsets] = msptdpcref_beat_detector(sig,fs, options)
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

% version: no optimisation

%% Options for running algorithm
if nargin<3
    options = struct;
end
options = setup_default_options(options);

%% Window signal
% split into overlapping windows
no_samps_in_win = options.win_len*fs;
if length(sig) <= no_samps_in_win
    win_starts = 1;
    win_ends = length(sig);
else
    win_offset = round(no_samps_in_win*(1-options.win_overlap));
    win_starts = 1:win_offset:length(sig)-no_samps_in_win;
    win_ends = win_starts + no_samps_in_win;
    if win_ends(end) < length(sig)
        win_starts(end+1) = length(sig) - no_samps_in_win;
        win_ends(end+1) = length(sig);
    end % this ensures that the windows include the entire signal duration
end

%% Downsample signal

% Set up downsampling if the sampling frequency is particularly high
if options.do_ds
    min_fs = options.ds_freq;
    if fs > min_fs
        ds_factor = floor(fs/min_fs);
        ds_fs = fs/floor(fs/min_fs);
    else
        options.do_ds = 0;
    end
end

%% detect peaks and onsets in each window

peaks = []; onsets = [];
for win_no = 1 : length(win_starts)

    % - extract this window's data
    win_sig = sig(win_starts(win_no):win_ends(win_no));
    
    % - downsample signal
    if options.do_ds
        rel_sig = downsample(win_sig, ds_factor);
        rel_fs = ds_fs;
    else
        rel_sig = win_sig;
        rel_fs = fs;
    end
    
    % - detect peaks and onsets
    [p,t] = detect_peaks_and_onsets_using_msptd(rel_sig, rel_fs, options);
    
    % - resample peaks
    if options.do_ds
        p = p*ds_factor;
        t = t*ds_factor;
    end
    
    % - correct peak indices by finding highest point within tolerance either side of detected peaks
    tol_durn = 0.05;
    if rel_fs < 10
        tol_durn = 0.2;
    elseif rel_fs < 20
        tol_durn = 0.1;
    end
    tol = ceil(fs*tol_durn);
    for pk_no = 1 : length(p)
        [~, temp] = max(win_sig( (p(pk_no) - tol) : (p(pk_no) + tol) ));
        p(pk_no) = p(pk_no) - tol + temp - 1;
    end
    
    % - correct onset indices by finding highest point within tolerance either side of detected onsets
    for onset_no = 1 : length(t)
        [~, temp] = min(win_sig( (t(onset_no) - tol) : (t(onset_no) + tol) ));
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

% % (i) to account for downsampling
% if options.do_ds
%     % - peaks
%     peaks = peaks*ds_factor;
%     % find highest point within tolerance either side of detected peaks
%     for pk_no = 1 : length(peaks)
%         curr_peak = peaks(pk_no);
%         tol_start = curr_peak - ds_factor;
%         tol_end = curr_peak + ds_factor;
%         [~, temp] = max(sig(tol_start:tol_end));
%         peaks(pk_no) = curr_peak - ds_factor + temp - 1;
%         clear temp curr_peak tol_start tol_end
%     end
%     % - onsets
%     onsets = onsets*ds_factor;
%     % find highest point within tolerance either side of detected onsets
%     for onset_no = 1 : length(onsets)
%         curr_onset = onsets(onset_no);
%         tol_start = curr_onset - ds_factor;
%         tol_end = curr_onset + ds_factor;
%         [~, temp] = min(sig(tol_start:tol_end));
%         onsets(onset_no) = curr_onset - ds_factor + temp - 1;
%         clear temp curr_onset tol_start tol_end
%     end
% end

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

function [p,t] = detect_peaks_and_onsets_using_msptd(x,fs, options)

% Setup
N = length(x); % length of signal
L = ceil(N/2)-1; % max window length

% Step 0: (added by PC) don't calculate scales outside the range of plausible HRs
plaus_hr_hz = options.plaus_hr_bpm./60; % in Hz
init_scales = 1:L;
durn_signal = length(x)/fs;
init_scales_fs = (L./init_scales)/durn_signal;
%init_scales_inc_log = init_scales_fs >= plaus_hr_hz(1) & init_scales_fs <= plaus_hr_hz(2);
if options.use_reduced_lms_scales
    init_scales_inc_log = init_scales_fs >= plaus_hr_hz(1);
else
    init_scales_inc_log = true(size(init_scales_fs));
end
max_scale = find(init_scales_inc_log,1,'last');

% Step 1: calculate local maxima and local minima scalograms

% - detrend
x = detrend(x); % this removes the best-fit straight line

% - initialise LMS matrices
% m_max = false(L,N);
% m_min = false(L,N);
if options.find_pks
    m_max = false(max_scale,N);
end
if options.find_trs
    m_min = false(max_scale,N);
end

% - populate LMS matrices
if options.optimisation
    options.lms_calc_method = 6;
end

switch options.lms_calc_method
    case 6

        % decide whether to use pks or trs
        if options.find_pks
            pk_log = true;
        else
            pk_log = false;
        end

        % Define anonymous function
        fun = @(k) calc_inv_gamma(k, max_scale, x, pk_log);
        
        % Find the value of k that results in a minimum value of the inverse of gamma
        fminbnd_options = optimset('MaxFunEvals',10,'TolX', 1); %, 'Display','iter');
        %max_scale = fminbnd(fun, 1, max_scale, fminbnd_options);
        % - suppress output
        [~, max_scale] = evalc('fminbnd(fun, 1, max_scale, fminbnd_options);');
        max_scale = ceil(max_scale);

        [m_max, m_min] = find_lms_using_msptd_approach(max_scale, x, options);

    case 5
        % Create indices for each scale
        indices = repmat((1:N), max_scale, 1);
        
        % Define indices for differences
        diff_indices1 = indices + (1:max_scale)';
        diff_indices2 = indices - (1:max_scale)';
        
        % Ensure the indices are within bounds
        bad_els = diff_indices1 > N | diff_indices2 < 1;
        diff_indices1(bad_els) = indices(bad_els);
        diff_indices2(bad_els) = indices(bad_els);
        
        % Calculate differences using vectorized operations
        %diff_matrix1 = x(indices) - x(diff_indices1);
        %diff_matrix2 = x(indices) - x(diff_indices2);
        S1 = x(indices);
        S2 = x(diff_indices2);
        S3 = x(diff_indices1);
        diff_matrix1_g = S1 > S3;
        diff_matrix2_g = S1 > S2;
        diff_matrix1_l = S1 < S3;
        diff_matrix2_l = S1 < S2;
        
        % populate
        m_max = diff_matrix1_g & diff_matrix2_g;
        m_min = diff_matrix1_l & diff_matrix2_l;
        
    case 4
        % Create indices for each scale
        tic
        indices = repmat((1:N), max_scale, 1);
        toc

        % Define indices for differences
        tic
        diff_indices1 = indices + (1:max_scale)';
        diff_indices2 = indices - (1:max_scale)';
        toc

        % Ensure the indices are within bounds
        tic
        diff_indices1(diff_indices1 > N) = indices(diff_indices1 > N);
        diff_indices2(diff_indices2 < 1) = indices(diff_indices2 < 1);
        toc

        % Calculate differences using vectorized operations
        %diff_matrix1 = x(indices) - x(diff_indices1);
        %diff_matrix2 = x(indices) - x(diff_indices2);
        tic
        diff_matrix1_g = x(indices) > x(diff_indices1);
        diff_matrix2_g = x(indices) > x(diff_indices2);
        diff_matrix1_l = x(indices) < x(diff_indices1);
        diff_matrix2_l = x(indices) < x(diff_indices2);
        toc

        % populate
        tic
        m_max = diff_matrix1_g & diff_matrix2_g;
        m_min = diff_matrix1_l & diff_matrix2_l;
        toc

    case 3
        diff_matrix1 = zeros(size(m_max));
        diff_matrix2 = zeros(size(m_max));
        for k = 1:max_scale % scalogram scales
            S1_start_el = k+1;
            S1_end_el = N-k-1;
            S1_els = S1_start_el:S1_end_el;
            diff_matrix1(k,S1_els) = x(S1_els) - x(S1_els-k);
            diff_matrix2(k,S1_els) = x(S1_els) - x(S1_els+k);
        end
        % populate for each scale, k, (i.e. row) in turn
        m_max = diff_matrix1>0 & diff_matrix2>0;
        m_min = diff_matrix1<0 & diff_matrix2<0;

    case 2
        % populate for each scale, k, (i.e. row) in turn
        for k = 1:max_scale % scalogram scales
            S1_start_el = k+1;
            S1_end_el = N-k-1;
            S1_els = S1_start_el:S1_end_el;
            S1 = x(S1_els);
            S2 = x(S1_els-k);
            S3 = x(S1_els+k);
            m_max(k,S1_els) = S1>S2 & S1>S3;
            m_min(k,S1_els) = S1<S2 & S1<S3;
        end
    case 1

        [m_max, m_min] = find_lms_using_msptd_approach(max_scale, x, options);

end

% Step 2: find the scale with the most local maxima (or local minima)
% - row-wise summation (i.e. sum each row)
if options.find_pks
    gamma_max = sum(m_max,2); % the ",2" option makes it row-wise
end
if options.find_trs
    gamma_min = sum(m_min,2); % the ",2" option makes it row-wise
end
% - find scale with the most local maxima (or local minima)
if options.find_pks
    [~, lambda_max] = max(gamma_max);
end
if options.find_trs
    [~, lambda_min] = max(gamma_min);
end

% Step 3: Use lambda to remove all elements of m for which k>lambda
first_scale_to_include = find(init_scales_inc_log,1);
if options.find_pks
    m_max = m_max(first_scale_to_include:lambda_max,:);
end
if options.find_trs
    m_min = m_min(first_scale_to_include:lambda_min,:);
end

% Step 4: Find peaks
% - column-wise summation
if options.find_pks
    m_max_sum = sum(~m_max);
    p = find(~m_max_sum); p = p(:);
else
    p = [];
end
if options.find_trs
    m_min_sum = sum(~m_min);
    t = find(~m_min_sum); t = t(:);
else
    t = [];
end
end

function inv_gamma = calc_inv_gamma(k, max_scale, x, pk_log)

k = round(k);

m = false(1,max_scale);
N = length(x);
for i = (k+2):(N-k+1)
    % maxima
    if pk_log
        if x(i-1) > x(i-k-1) && x(i-1) > x(i+k-1)
            m(1,i-1) = true;
        end
    else
        % minima
        if x(i-1) < x(i-k-1) && x(i-1) < x(i+k-1)
            m(1,i-1) = true;
        end
    end
end

gamma = sum(m);

inv_gamma = 1./gamma;

end

function [m_max, m_min] = find_lms_using_msptd_approach(max_scale, x, options)

N = length(x);

m_max = false(max_scale,N);
m_min = false(max_scale,N);
        
if options.find_pks && options.find_trs
    
    m_max = find_m_max(x, N, max_scale, m_max);

    m_min = find_m_min(x, N, max_scale, m_min);

elseif options.find_pks
    
    m_max = find_m_max(x, N, max_scale, m_max);

elseif options.find_trs

    m_min = find_m_min(x, N, max_scale, m_min);

end

end

function m_max = find_m_max(x, N, max_scale, m_max)

% populate for each scale, k, (i.e. row) in turn
for k = 1:max_scale % scalogram scales
    % populate for each column
    for i = (k+2):(N-k+1)
        % maxima
        if x(i-1) > x(i-k-1) && x(i-1) > x(i+k-1)
            m_max(k,i-1) = true;
        end
    end
end

end

function m_min = find_m_min(x, N, max_scale, m_min)

% populate for each scale, k, (i.e. row) in turn
for k = 1:max_scale % scalogram scales
    % populate for each column
    for i = (k+2):(N-k+1)
        % minima
        if x(i-1) < x(i-k-1) && x(i-1) < x(i+k-1)
            m_min(k,i-1) = true;
        end
    end
end

end

function options = setup_default_options(options)
% Setups up the options for the analysis, using the values provided by the 
% user, and default values for options not provided by the user.

%% Specify all the options
option_names = {'find_pks'; 'find_trs'; 'lms_calc_method'; 'use_reduced_lms_scales'; 'do_ds'; 'ds_freq'; 'win_len'; 'win_overlap'; 'optimisation'; 'plaus_hr_bpm'};

%% Specify the setting for each option

% cycle through each option in turn
for s = 1 : length(option_names)
    % - check to see whether a setting has been provided
    if ~sum(strcmp(fieldnames(options),option_names{s}))
        
        % if not, then use default option
        switch option_names{s}
            case 'find_pks'
                default_setting = true;
            case 'find_trs'
                default_setting = true;
            case 'lms_calc_method'
                default_setting = 1;
            case 'use_reduced_lms_scales'
                default_setting = false;
            case 'do_ds'
                default_setting = false;
            case 'ds_freq'
                default_setting = nan;
            case 'win_len'
                default_setting = 8;
            case 'win_overlap'
                default_setting = 0.2;
            case 'optimisation'
                default_setting = false;
            case 'plaus_hr_bpm'
                default_setting = [30, 200];
        end
        % store this setting in the uParams.analysis structure:
        eval(['options.' option_names{s} ' = default_setting;']);
        clear default_setting
    end
end

end


