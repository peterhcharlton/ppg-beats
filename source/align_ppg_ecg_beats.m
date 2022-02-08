function [ecg_beats_a_inds, ecg_exc_a_log] = align_ppg_ecg_beats(ppg_beats_inds, ecg_beats_inds, ppg_fs, ecg_fs, options, ecg_exc_log)
% ALIGN_PPG_ECG_BEATS  aligns PPG and ECG beats.
%   ALIGN_PPG_ECG_BEATS aligns beats detected in an electrocardiogram (ECG)
%   signal with those detected in a simultaneous photoplethysmogram (PPG) signal.
%   
%   # Inputs
%   
%   * ppg_beats_inds  - indices of PPG beats
%   * ecg_beats_inds  - indices of ECG beats
%   * ppg_fs          - sampling frequency of PPG in Hz
%   * ecg_fs          - sampling frequency of ECG in Hz
%   * options         - (optional) A structure containing the settings for the analysis:
%    - options.max_lag     - max permissible lag between ECG and PPG in secs
%    - options.lag_int     - increment between tested lags
%    - options.tol_window  - acceptable tolerance between ECG and PPG beats
%   * ecg_exc_log     - (optional) A logical indicating whether or not each ECG sample will be excluded from the analysis
%   
%   # Outputs
%   * ecg_beats_a_inds : the aligned indices of ECG beats
%   * ecg_exc_a_log : the aligned ECG exclusion logical
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Author
%   Peter H. Charlton, University of Cambridge, February 2022.
%   
%   # Version
%   0.1, and is still in development.
%   
%   # Licence
%      This file is part of PPG-beats.
%      PPG-beats is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%      PPG-beats is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%      You should have received a copy of the GNU General Public License along with PPG-beats. If not, see <https://www.gnu.org/licenses/>.

% in the future this could be adjusted to include exclusion of low quality ECG beats with ecg_exc_log.

%% Setup 

% - setup default options
options = setup_default_options(options);

% - skip if there are no detected PPG beats
if isempty(ppg_beats_inds)
    ecg_beats_a_inds = ecg_beats_inds;
    ecg_exc_a_log = ecg_exc_log;
    return
end

%% Create time vectors of beat indices
ecg_beats_t = [ecg_beats_inds-1]./ecg_fs; % (starting at time of 0)
ppg_beats_t = [ppg_beats_inds-1]./ppg_fs;

do_windows = 0;
if ~do_windows
    
    %% Assess how well beat detections agree for different time lags between the signals
    
    % - Calculate time differences between each pair of PPG and ECG beats
    diff_matrix = repmat(ppg_beats_t, [1, length(ecg_beats_t)]) - ecg_beats_t';
    
    % - define the lags to be tried
    lag_ts = (-1*options.max_lag):options.lag_int:options.max_lag;
    
    % - calculate performance for a range of lags
    curr_perf = nan(size(lag_ts));
    for lag_t_el = 1:length(lag_ts)
        
        % - identify the current lag
        lag_t = lag_ts(lag_t_el);
        
        % - Find the minimum difference between each ECG beat and its nearest PPG beat when using this lag
        min_abs_diff = min(abs(diff_matrix-lag_t));
        min_abs_diff = min_abs_diff(:);
        
        % - Determine whether or not each ECG beat is within an acceptable tolerance of a PPG beat when using this lag
        beat_correct_log = min_abs_diff<options.tol_window;
        
        % - Calculate performance as the number of ECG beats within an acceptable tolerance
        curr_perf(lag_t_el) = sum(beat_correct_log);
    end
    clear min_abs_diff beat_correct_log lag_t lag_t_el
    
    %% identify the lag which gives the optimal performance
    temp = max(curr_perf);
    rel_lag_t_els = find(curr_perf == temp);        % because more than one lag can result in this performance
    [~, temp2] = min(abs(lag_ts(rel_lag_t_els)));   % identify the index of the minimum absolute lag which results in this optimal performance
    rel_lag_t = lag_ts(rel_lag_t_els(temp2));       % identify the lag which results in the optimal performance
    rel_lag_ecg_samps = round(rel_lag_t*ecg_fs);    % express the lag in terms of the number of ECG samples which corresponds to this lag
    clear temp2 temp curr_perf rel_lag_t_els lag_ts
    
    %% calculate the aligned indices of ECG beats
    ecg_beats_a_inds = ecg_beats_inds+rel_lag_ecg_samps;
    
    %% calculate time-shifted ecg exclusion log
    if exist('ecg_exc_log', 'var')
        if rel_lag_ecg_samps > 0
            ecg_exc_a_log = [true(rel_lag_ecg_samps,1); ecg_exc_log(1:end-rel_lag_ecg_samps)];
        elseif rel_lag_ecg_samps < 0
            ecg_exc_a_log = [ecg_exc_log(abs(rel_lag_ecg_samps)+1:end); true(abs(rel_lag_ecg_samps),1)];
        else
            ecg_exc_a_log = ecg_exc_log;
        end
    else
        ecg_exc_a_log = [];
    end
    
else
    
    %% Split into windows
    win_durn = 120; % in secs
    win_durn_samps = win_durn*ecg_fs;
    win_starts = 1:win_durn_samps:(length(ecg_exc_log)-win_durn_samps);
    
    ecg_beats_a_inds = [];
    if exist('ecg_exc_log', 'var')
        ecg_exc_a_log = true(size(ecg_exc_log));
    end
    
    for win_no = 1 : length(win_starts)
        curr_win_start = win_starts(win_no);
        curr_win_end = curr_win_start+win_durn_samps-1;
        curr_rel_ecg_beats = ecg_beats_inds>=curr_win_start & ecg_beats_inds<curr_win_end;
        curr_ecg_beats_t = ecg_beats_t(curr_rel_ecg_beats);
        curr_rel_ppg_beats = ppg_beats_inds>=curr_win_start & ppg_beats_inds<curr_win_end;
        curr_ppg_beats_t = ppg_beats_t(curr_rel_ppg_beats);
        
        no_ecg_samps = find_lag_between_ecg_and_ppg(curr_ppg_beats_t, curr_ecg_beats_t, ecg_fs, options);
        %fprintf('\n - no samps: %d', no_ecg_samps);
        
        %% calculate the aligned indices of ECG beats
        curr_ecg_beats_a_inds = ecg_beats_inds(curr_rel_ecg_beats)+no_ecg_samps;
        if ~isempty(ecg_beats_a_inds)
            curr_ecg_beats_a_inds = curr_ecg_beats_a_inds(curr_ecg_beats_a_inds > ecg_beats_a_inds(end));
        end
        ecg_beats_a_inds = [ecg_beats_a_inds; curr_ecg_beats_a_inds];
        
        %% calculate time-shifted ecg exclusion log
        if exist('ecg_exc_log', 'var')
            start_el = curr_win_start + no_ecg_samps;
            end_el = curr_win_end + no_ecg_samps;
            if start_el < 1
                no_to_exclude_at_start = 1-start_el;
            else
                no_to_exclude_at_start = 0;
            end
            if end_el > length(ecg_exc_log)
                no_to_exclude_at_end = end_el - length(ecg_exc_log);
            else
                no_to_exclude_at_end = 0;
            end
            ecg_exc_a_log(start_el+no_to_exclude_at_start:end_el-no_to_exclude_at_end) = ecg_exc_log(curr_win_start+no_to_exclude_at_start:curr_win_end-no_to_exclude_at_end);
        end
        
    end
    
end

% To check:
% plot(ppg_beats_t, zeros(size(ppg_beats_t)), '*k'), hold on, plot(ecg_beats_t, zeros(size(ecg_beats_t)), 'or'), plot(ppg_beats_t, ones(size(ppg_beats_t)), '*k'), hold on, plot(ecg_beats_t+rel_lag_t, ones(size(ecg_beats_t)), 'or'), ylim([-1, 2])

end

function options = setup_default_options(options)
% Setups up the options for the analysis, using the values provided by the 
% user, and default values for options not provided by the user.

%% Specify all the options
option_names = {'max_lag'; 'lag_int'; 'tol_window'};

%% Specify the setting for each option

% cycle through each option in turn
for s = 1 : length(option_names)
    % - check to see whether a setting has been provided. If not, then use the default
    if ~sum(strcmp(fieldnames(options),option_names{s}))
        switch option_names{s}
            case 'max_lag'
                default_setting = 10;  % max permissible lag between ECG and PPG in secs
            case 'lag_int'
                default_setting = 0.02;
            case 'tol_window'
                default_setting = 0.2; % in secs
        end
        % store this setting in the options structure:
        eval(['options.' option_names{s} ' = default_setting;']);
        clear default_setting
    end
end
clear s

end

function no_ecg_samps = find_lag_between_ecg_and_ppg(ppg_beats_t, ecg_beats_t, ecg_fs, options)


%% Assess how well beat detections agree for different time lags between the signals

% - Calculate time differences between each pair of PPG and ECG beats
diff_matrix = repmat(ppg_beats_t, [1, length(ecg_beats_t)]) - ecg_beats_t';

% - define the lags to be tried
lag_ts = (-1*options.max_lag):options.lag_int:options.max_lag;

% - calculate performance for a range of lags
curr_perf = nan(size(lag_ts));
for lag_t_el = 1:length(lag_ts)
    
    % - identify the current lag
    lag_t = lag_ts(lag_t_el);
    
    % - Find the minimum difference between each ECG beat and its nearest PPG beat when using this lag
    min_abs_diff = min(abs(diff_matrix-lag_t));
    min_abs_diff = min_abs_diff(:);
    
    % - Determine whether or not each ECG beat is within an acceptable tolerance of a PPG beat when using this lag
    beat_correct_log = min_abs_diff<options.tol_window;
    
    % - Calculate performance as the number of ECG beats within an acceptable tolerance
    curr_perf(lag_t_el) = sum(beat_correct_log);
end
clear min_abs_diff beat_correct_log lag_t lag_t_el

%% identify the lag which gives the optimal performance
temp = max(curr_perf);
rel_lag_t_els = find(curr_perf == temp);        % because more than one lag can result in this performance
[~, temp2] = min(abs(lag_ts(rel_lag_t_els)));   % identify the index of the minimum absolute lag which results in this optimal performance
rel_lag_t = lag_ts(rel_lag_t_els(temp2));       % identify the lag which results in the optimal performance
no_ecg_samps = round(rel_lag_t*ecg_fs);         % express the lag in terms of the number of ECG samples which corresponds to this lag
clear temp2 temp curr_perf rel_lag_t_els lag_ts

end