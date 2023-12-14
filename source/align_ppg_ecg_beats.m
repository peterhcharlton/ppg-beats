function [ecg_beats_a_inds, ecg_exc_a_log, lag_ecg_samps] = align_ppg_ecg_beats(ppg_beats_inds, ecg_beats_inds, ppg_fs, ecg_fs, options, ecg_exc_log)
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
%    - options.do_wins     - whether or not to use windowing to find the time lag
%   * ecg_exc_log     - (optional) A logical indicating whether or not each ECG sample will be excluded from the analysis
%   
%   # Outputs
%   * ecg_beats_a_inds : the aligned indices of ECG beats
%   * ecg_exc_a_log :    the aligned ECG exclusion logical
%   * lag_ecg_samps :    the number of samples by which the ECG is in front / behind the PPG
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Author
%   Peter H. Charlton, University of Cambridge, February 2022.
%   
%   # Version
%   1.1
%   
%   # License - GPL-3.0
%      Copyright (c) 2022 Peter H. Charlton
%      This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%      This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%      You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

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

%% Identify lag between PPG and ECG beats
lag_ecg_samps = identify_lag_between_ppg_and_ecg(ppg_beats_t, ecg_beats_t, ecg_fs, options);

%% calculate the aligned indices of ECG beats
ecg_beats_a_inds = ecg_beats_inds+lag_ecg_samps;
ecg_beats_a_inds(ecg_beats_a_inds<0) = []; % some are shifted to before time 0 due to the lag.

%% calculate time-shifted ecg exclusion log
if exist('ecg_exc_log', 'var')
    if unique(lag_ecg_samps) == 0
        ecg_exc_a_log = ecg_exc_log;
    elseif length(lag_ecg_samps) == 1 % i.e. for non-variable lag
        if lag_ecg_samps > 0
            ecg_exc_a_log = [true(lag_ecg_samps,1); ecg_exc_log(1:end-lag_ecg_samps)];
        elseif lag_ecg_samps < 0
            ecg_exc_a_log = [ecg_exc_log(abs(lag_ecg_samps)+1:end); true(abs(lag_ecg_samps),1)];
        end
    else % i.e. for variable lag
        % shift exclusion log
        old_sig_inds = 1:length(ecg_exc_log);
        new_sig_lags = round(interp1(ecg_beats_inds,lag_ecg_samps,old_sig_inds));
        new_sig_inds = old_sig_inds+new_sig_lags;
        new_sig_inds(new_sig_inds<1) = nan;
        ecg_exc_a_log = nan(length(ecg_exc_log),1); % setup variable
        ecg_exc_a_log(new_sig_inds(~isnan(new_sig_inds))) = ecg_exc_log(~isnan(new_sig_inds)); % insert old exc log at lagged points in new one
        ecg_exc_a_log = ecg_exc_a_log(1:length(ecg_exc_log)); % truncate any overlength indices
        ecg_exc_a_log(isnan(ecg_exc_a_log)) = 1;

        % exclude any beats which are not monotonically increasing in time
        curr_max = 0;
        non_monotonic_beats = false(length(ecg_beats_a_inds),1);
        for beat_no = 1 : length(ecg_beats_a_inds)
            if ecg_beats_a_inds(beat_no)>curr_max
                curr_max = ecg_beats_a_inds(beat_no);
            else
                non_monotonic_beats(beat_no) = true;
            end
        end
        if sum(non_monotonic_beats)
            non_monotonic_beat_inds = ecg_beats_a_inds(non_monotonic_beats);
            % mark these for exclusion
            ecg_exc_a_log(non_monotonic_beat_inds) = true;
            % exclude these beats
            ecg_beats_a_inds(non_monotonic_beats) = [];
        end
        
    end
else
    ecg_exc_a_log = [];
end

% To check:
% plot(ppg_beats_t, zeros(size(ppg_beats_t)), '*k'), hold on, plot(ecg_beats_t, zeros(size(ecg_beats_t)), 'or'), plot(ppg_beats_t, ones(size(ppg_beats_t)), '*k'), hold on, plot(ecg_beats_t+rel_lag_t, ones(size(ecg_beats_t)), 'or'), ylim([-1, 2])

end

function lag_ecg_samps = identify_lag_between_ppg_and_ecg(ppg_beats_t, ecg_beats_t, ecg_fs, options)

%% Assess how well beat detections agree for different time lags between the signals

% - Calculate time differences between each pair of PPG and ECG beats
diff_matrix = repmat(ppg_beats_t, [1, length(ecg_beats_t)]) - ecg_beats_t';

% - define the lags to be tried
lag_ts = (-1*options.max_lag):options.lag_int:options.max_lag;

% - calculate performance for a range of lags
if options.do_wins
    win_durn = 60; % in secs
    win_starts = 0:win_durn:ecg_beats_t(end);
    win_starts(end+1) = ecg_beats_t(end)-win_durn;
    win_ends = win_starts + win_durn;
end
if options.do_variable_lag
    all_diffs = nan(length(ecg_beats_t), length(lag_ts));
    lag_win_durn = 300;
    lag_win_starts = 1:lag_win_durn:length(ecg_beats_t)-lag_win_durn;
    if isempty(lag_win_starts), lag_win_starts = 1; end  % CHANGE added by PC on 3-Oct-2023 for short signal duration
    lag_win_ends = lag_win_starts + lag_win_durn - 1;
    lag_win_ends(end) = length(ecg_beats_t);
    curr_perf = nan(length(lag_ts), length(lag_win_starts));
else
    curr_perf = nan(size(lag_ts));
end

for lag_t_el = 1:length(lag_ts)
    
    % - identify the current lag
    lag_t = lag_ts(lag_t_el);
    
    % - Find the minimum difference between each ECG beat and its nearest PPG beat when using this lag
    if options.do_wins
        min_abs_diff = nan(length(ecg_beats_t),1);
        for win_no = 1 : length(win_starts)
            ppg_els = ppg_beats_t>=(win_starts(win_no)-options.max_lag) & ppg_beats_t<=(win_ends(win_no)+options.max_lag);
            ecg_els = ecg_beats_t>=win_starts(win_no) & ecg_beats_t<=win_ends(win_no);
            if ~sum(ppg_els) || ~sum(ecg_els)
                continue
            end
            min_abs_diff(ecg_els) = find_min_abs_diff(diff_matrix(ppg_els,ecg_els)-lag_t);
        end
        if options.do_variable_lag
            all_diffs(:,lag_t_el) = min_abs_diff;
        end
    else
        min_abs_diff = find_min_abs_diff(diff_matrix-lag_t);
    end
    
    % - Determine whether or not each ECG beat is within an acceptable tolerance of a PPG beat when using this lag
    beat_correct_log = min_abs_diff<options.tol_window;
    if options.do_variable_lag
        all_beat_correct_log = all_diffs<options.tol_window;
    end
    
    % - Calculate performance as the number of ECG beats within an acceptable tolerance
    if options.do_variable_lag
        for lag_win_no = 1 : length(lag_win_starts)
            curr_perf(lag_t_el,lag_win_no) = sum(beat_correct_log(lag_win_starts(lag_win_no):lag_win_ends(lag_win_no)));
        end
    else
        curr_perf(lag_t_el) = sum(beat_correct_log);
    end
end
clear min_abs_diff beat_correct_log lag_t lag_t_el

%% identify the lag which gives the optimal performance
if options.do_variable_lag
    lag_ecg_samps = nan(length(ecg_beats_t),1);
    for lag_win_no = 1 : length(lag_win_starts)
        rel_curr_perf = curr_perf(:,lag_win_no);
        win_lag_ecg_samps = find_lag_for_curr_perf(rel_curr_perf, lag_ts, ecg_fs);
        lag_ecg_samps(lag_win_starts(lag_win_no):lag_win_ends(lag_win_no)) = win_lag_ecg_samps;
    end
else
    lag_ecg_samps = find_lag_for_curr_perf(curr_perf, lag_ts, ecg_fs);    
end

end

function lag_ecg_samps = find_lag_for_curr_perf(curr_perf, lag_ts, ecg_fs)

temp = max(curr_perf);
rel_lag_t_els = find(curr_perf == temp);        % because more than one lag can result in this performance
[~, temp2] = min(abs(lag_ts(rel_lag_t_els)));   % identify the index of the minimum absolute lag which results in this optimal performance
rel_lag_t = lag_ts(rel_lag_t_els(temp2));       % identify the lag which results in the optimal performance
lag_ecg_samps = round(rel_lag_t*ecg_fs);        % express the lag in terms of the number of ECG samples which corresponds to this lag

end

function min_abs_diff = find_min_abs_diff(diff_matrix)

% find minimum in each column of 'diff_matrix'

min_abs_diff = min(abs(diff_matrix));
min_abs_diff = min_abs_diff(:);
    
end

function lag_ecg_samps = identify_lag_between_ppg_and_ecg_old(ppg_beats_t, ecg_beats_t, options);

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
lag_ecg_samps = round(rel_lag_t*ecg_fs);        % express the lag in terms of the number of ECG samples which corresponds to this lag

%clear temp2 temp curr_perf rel_lag_t_els lag_ts

end

function options = setup_default_options(options)
% Setups up the options for the analysis, using the values provided by the 
% user, and default values for options not provided by the user.

%% Specify all the options
option_names = {'max_lag'; 'lag_int'; 'tol_window'; 'do_wins'; 'do_variable_lag'};

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
            case 'do_wins'
                default_setting = true;
            case 'do_variable_lag'
                default_setting = true;
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