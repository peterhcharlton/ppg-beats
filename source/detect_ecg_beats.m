function [beats_inds, qual] = detect_ecg_beats(ecg, fs, options, no_signal_vector)
% DETECT_ECG_BEATS  detects beats in ECG.
%   DETECT_ECG_BEATS detects beats in an electrocardiogram (ECG) signal
%   using two beat detectors, and uses the agreement between beat detectors to
%   assess the quality of beat detections.
%   
%   # Inputs
%   
%   * ecg : a vector of ECG values
%   * fs : the sampling frequency of the ECG in Hz
%   * options : (optional) a structure of options which determine the settings used for the analysis:
%    - options.win_durn       - The duration of windows over which to perform quality assessment (in secs)
%    - options.win_overlap    - The overlap of windows (in secs)
%    - options.verbose        - A logical indicating whether (1) or not (0) to display text in the Command Window
%    - options.qrs_tol_window - The acceptable tolerance window around QRS spikes (in secs)
%
%   * no_signal_vector : (optional) a vector of the same length as the ECG vector, which indicates whether (1) or not (0) there was a signal at each sample
%   
%   # Outputs
%   
%   * beat_inds : the indices of detected beats which both beat detectors agree on
%   * qual : a logical indicating whether (1) or not (0) beats agreed between beat detectors for each ECG sample 
%   * temporary file - this script creates a temporary file in the current directory
%   
%   # Exemplary usage:
%
%       options.win_durn = 20;
%       [beat_inds, qual] = detect_ecg_beats(ecg, fs, options)
%   
%   # Requirements:
%   
%   Waveform Database Software Package (WFDB) for MATLAB and Octave, available from:
%   <https://doi.org/10.13026/6zcz-e163>
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Author
%   Peter H. Charlton, University of Cambridge, February 2022.
%
%   # Acknowledgment:
%   This script uses scripts from the PhysioNet Cardiovascular Signal Toolbox (under the GNU public licence):
%   - Scripts: ExportHRVparams.m, run_qrsdet_by_seg.m, jqrs.m
%   - Some of the parameters included in this script are taken from InitializeHRVparams.m
%   - Link: <https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox>
%   - Reference: Vest et al. 'An open source benchmarked toolbox for cardiovascular waveform and interval analysis'. Physiol Meas, 39, 2018. <https://doi.org/10.1088/1361-6579/aae021> 
%   
%   # Version
%   0.1, and is still in development.
%   
%   # Licence
%      This file is part of PPG-beats.
%      PPG-beats is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%      PPG-beats is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%      You should have received a copy of the GNU General Public License along with PPG-beats. If not, see <https://www.gnu.org/licenses/>.

%% Setup

% - ECG vector
ecg = ecg(:); % make into column vector

% - replace any nans with median value of ECG
ecg(isnan(ecg)) = median(ecg(~isnan(ecg)));

% - options
if ~exist('options', 'var')
    options = struct;
end
options = setup_options(options);

% - Display startup message
if options.start_up_message, display_startup_message(options); end

% - WFDB toolbox
check_WFDB_installation;

%% Detect beats using the GQRS detector

if options.verbose, fprintf('\n - Detecting beats using the gqrs beat detector: '), end

% - write this signal to a WFDB file
%   - identify apparent resolution
temp_diffs = diff(sort(ecg));
temp_sig = round(ecg./(min(temp_diffs(temp_diffs>0))));
%   - create time vector
tm = (0:(length(ecg)-1))/fs;
%   - write to file
wrsamp(tm(:),temp_sig,'temp',fs);

% - run gqrs to identify beats, and save to a temporary qrs annotation file
gqrs('temp');

% - load qrs annotations
[gqrs_beat_inds]=rdann('temp', 'qrs');
if options.verbose, fprintf('%d beats detected', length(gqrs_beat_inds)), end

% - display message about temporary files
if options.verbose
    fprintf('\n   - temporary files created at: %s%stemp', cd, filesep)
end

%% JQRS detector

if options.verbose, fprintf('\n - Detecting beats using the jqrs beat detector: '), end

% - Set parameters
HRVparams.PeakDetect.REF_PERIOD = 0.250; 
HRVparams.PeakDetect.THRES = .6; 
HRVparams.Fs = fs;
HRVparams.PeakDetect.fid_vec = [];
HRVparams.PeakDetect.SIGN_FORCE = [];
HRVparams.PeakDetect.debug = 0;
HRVparams.PeakDetect.windows = 15;
HRVparams.PeakDetect.ecgType = 'MECG';

jqrs_beat_inds = run_qrsdet_by_seg(ecg,HRVparams);
jqrs_beat_inds = jqrs_beat_inds(:);
if options.verbose, fprintf('%d beats detected', length(jqrs_beat_inds)), end

%% Identify 'correct' beats

if options.verbose, fprintf('\n - Identifying correct beat detections: '), end

if isempty(jqrs_beat_inds) || isempty(gqrs_beat_inds)
    beat_correct_log = false(size(jqrs_beat_inds));
else
    % - Calculate difference matrix using gqrs as the reference
    diff_matrix = repmat(jqrs_beat_inds, [1, length(gqrs_beat_inds)]) - gqrs_beat_inds';
    % - Find minimum differences
    min_abs_diff = min(abs(diff_matrix), [], 2);
    % - Identify correctly identified beats
    beat_correct_log = min_abs_diff<(options.qrs_tol_window*fs);
end

% identify beats as those which were agreed on by both beat detectors
beats_inds = jqrs_beat_inds(beat_correct_log);

% - Remove any repeated beat detections
beats_inds = sort(beats_inds);
if ~isempty(beats_inds)
    repeated_beats = [0; diff(beats_inds)< round(options.qrs_tol_window*fs)];
    beats_inds = beats_inds(~repeated_beats);
    clear repeated_beats
end
if options.verbose, fprintf('%d out of %d beats detected by gqrs deemed to be correct', length(beats_inds), length(gqrs_beat_inds)), end

%% Assess quality of beat detections in windows

if options.verbose, fprintf('\n - Assessing quality of beat detections in windows: '), end

% - Identify window start times
win_starts = 0:(options.win_durn-options.win_overlap):(((length(ecg)-1)/fs)-options.win_durn);
% - Create time vector
t = 0:(1/fs):((length(ecg)-1)/fs);
% - Identify times of disagreement between beat detectors
t_disagree_beats = t(jqrs_beat_inds(~beat_correct_log));
% - Determine quality of ECG in each window (i.e. if there was an incorrect beat detection, then deem it to be of low quality)
qual = false(length(ecg),1);
for win_no = 1 : length(win_starts)
    
    % - identify elements corresponding to this window
    curr_start = win_starts(win_no);
    curr_end = curr_start+options.win_durn;
    rel_els = t>= curr_start & t<= curr_end;
    
    % - skip this window if any of it didn't contain a signal
    if exist('no_signal_vector', 'var') && sum(no_signal_vector(rel_els))
        qual(rel_els) = true;
        continue
    end
    
    %  - deem this window to be of high quality if there were beats detected, and there weren't any disagreements between the beat detectors
    curr_win_disagree_inds = t_disagree_beats >= curr_start & t_disagree_beats <= curr_end;
    curr_win_beat_inds = t(jqrs_beat_inds) >= curr_start & t(jqrs_beat_inds) <= curr_end;
    if sum(curr_win_disagree_inds) == 0 && sum(curr_win_beat_inds)
        qual(rel_els) = true;
        continue
    end
    
    clear curr_start curr_end rel_els curr_win_disagree_inds
end

if options.verbose, fprintf('%.1f%% of ECG signal deemed to be of high quality', 100*mean(qual)), end

end

function options = setup_options(options)

specified_options = fieldnames(options);
possible_options = {'win_durn', 'win_overlap', 'verbose', 'qrs_tol_window', 'start_up_message'};

for option_no = 1 : length(possible_options)
    curr_option = possible_options{option_no};
    
    % skip if this option has already been specified
    if sum(strcmp(specified_options, curr_option))
        continue
    else
        % identify default value for this option
       switch curr_option
           case 'win_durn'
               default_val = 20; % default of 20 seconds
           case 'win_overlap'
               default_val = 5; % default of 5 seconds
           case 'verbose'
               default_val = true;
           case 'qrs_tol_window'
               default_val = 0.2;  % in secs
           case 'start_up_message'
               default_val = 1;    % i.e. display the start-up message
       end
       % store the default value for this option
        eval(['options.' curr_option ' = default_val;']);
    end
    
end

end

function check_WFDB_installation

% Check whether the 'wrsamp' function is available
if ~exist('wrsamp.m', 'file')
    error('This code requires the Waveform Database Software Package (WFDB) for MATLAB and Octave, available from: https://doi.org/10.13026/6zcz-e163');
end

end

function display_startup_message(options)
% Displays a starup message, including details of the licence

licence_details = ['\n\n DETECT_ECG_BEATS  Copyright (C) 2022  Peter H. Charlton',...
    '\n This program comes with ABSOLUTELY NO WARRANTY', ... 
    '\n and is available under the GNU public license.', ...
    '\n For details see the accompanying LICENSE.txt file.\n\n'];

if options.verbose
    fprintf(licence_details)
end

end