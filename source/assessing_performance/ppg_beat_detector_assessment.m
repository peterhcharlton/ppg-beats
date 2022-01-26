function ppg_beat_detector_assessment(dataset, options)
% PPG_BEAT_DETECTOR_ASSESSMENT  Assesses the performance of photoplethysmogram 
% (PPG) beat detectors. Performance can be assessed across a range of datasets.
%   
%  Inputs:
%
%     dataset - The name of the dataset with which to assess performance, e.g.
%                   'ppg_dalia_working'
%                   'ppg_dalia_walking'
%                   'ppg_dalia_sitting'
%                   'capnobase'
%                   'mimic_af'
%                   'mimic_non_af'
%                   'wesad_baseline'
%               You must have downloaded and collated the dataset in order
%               to use it. Instructions on how to download and collate datasets
%               are provided here:
%                   https://peterhcharlton.github.io/info/datasets
%
%     options - (optional) a structure of options which determine the settings used for the analysis:
%                   options.dataset_file    - a string containing the path of the Matlab file containing the collated dataset.
%                   options.beat_detectors  - a cell array containing strings corresponding to the names of the beat detectors to be used.
%                   options.do_downsample   - (1 or 0, default value of 0) A logical indicating whether or not to downsample the PPG signals prior to analysis.
%                   options.downsample_freq - the frequency at which to downsample PPG signals (Hz). Only used if downsampling is enabled by options.do_downsample.
%                   options.redo_analysis   - A logical indicating whether or not to overwrite existing analysis files in order to redo the analysis
%
%  Outputs:
%
%
%  Exemplary usage:
%
%       ppg_beat_detection_assessment('mimic_non_af')      % assesses the performance of beat detectors on the 'mimic_non_af' dataset using default analysis options
%
%  Further Information:
%       Please see the accompanying manual at: ...
%
%  Licence: 
%       Available under the GNU public license - please see the accompanying
%       file named "LICENSE"
%
% Author: Peter H. Charlton, University of Cambridge, July 2021.
% Version: 0.1, and is still in development.

%% Setup

% Initialise dataset and options if not specified:
if nargin == 0, dataset = ''; end
if nargin < 2, options = struct; end

% Setup universal parameters
uParams = setup_up(dataset, options);

uParams.analysis.redo_analysis = 0;

%% Detect beats in PPG signals
detect_beats_in_ppg_signals(uParams);

%% Assess quality of PPG signals
assess_quality_of_ppg_signals(uParams);

%% Detect beats in ECG signals
uParams.analysis.redo_analysis = 1;
detect_beats_in_ecg_signals(uParams);

%% Time align PPG beats
uParams.analysis.redo_analysis = 1;
time_align_ppg_beats(uParams);

% %% Assess performance of quality assessment tools
% assess_quality_assessment_tools(uParams);

%% Assess performance of PPG beat detectors
uParams.analysis.redo_analysis = 1;
assess_ppg_beat_detector_performance(uParams);

%% Calculate stats and create tables
uParams.analysis.redo_analysis = 1;
create_stats_and_tables(uParams);

end

function extra_stuff_to_incorporate

do_acc = 0;
if do_acc
    for beat_detector_no = 1 : length(uParams.analysis.beat_detectors)
        [acc_wins.mad{beat_detector_no}, acc_wins.no_ecg{beat_detector_no}, acc_wins.no_ppg{beat_detector_no}, acc_wins.no_correct{beat_detector_no}] = deal(nan(100000,1));
    end
end
[qual_sens, qual_spec, qual_bal_acc, sens, ppv, f1_score, no_beats, acc_hr, durn_hr, prop_hr, bias_hr, loa_hr] = deal(nan(length(uParams.analysis.beat_detectors), length(beats_t)));
for beat_detector_no = 1 : length(uParams.analysis.beat_detectors)
    
    acc_counter = 0;
    rel_lag_t = nan(length(beats_t),1);
    for subj_no = 1 : length(beats_t)
        
        %fprintf('\n   - Beat detector %s in subj %d:', up.analysis.beat_detectors{beat_detector_no}, subj_no);
                
        % --- quality assessment ---
        ppg_identified = beats_t(subj_no).ppg{beat_detector_no};
        diff_matrix = repmat(ppg_identified-rel_lag_t(subj_no), [1, length(beats_t(subj_no).ecg)]) - beats_t(subj_no).ecg';
        % - Find minimum differences
        min_abs_diff = min(abs(diff_matrix'));
        min_abs_diff = min_abs_diff(:);
        % - Identify correct beats
        beat_correct_log = min_abs_diff<uParams.analysis.tol_window;
        % - PPG quality values
        ppg_qual = beats_t(subj_no).ppg_qual{beat_detector_no};
        % - performance statistics
        qual_spec(beat_detector_no, subj_no) = 100*sum(~beat_correct_log & ~ppg_qual)/sum(~beat_correct_log);
        qual_sens(beat_detector_no, subj_no) = 100*sum(beat_correct_log & ppg_qual)/sum(beat_correct_log);
        qual_bal_acc(beat_detector_no, subj_no) = 0.5*(qual_sens(beat_detector_no, subj_no)+qual_spec(beat_detector_no, subj_no));
        
        clear diff_matrix min_abs_diff median_diff rel_ecg_els diff_matrix beat_correct_log ppg_identified ppg_qual
        
        
        % --- beat detection ---
        
        % - Calculate time differences between each pair of PPG and ECG beats
        hq_ppg_beats = beats_t(subj_no).ppg{beat_detector_no}(beats_t(subj_no).ppg_qual{beat_detector_no});
        diff_matrix = repmat(hq_ppg_beats-rel_lag_t(subj_no), [1, length(beats_t(subj_no).ecg)]) - beats_t(subj_no).ecg';
        % - Find minimum differences
        min_abs_diff = min(abs(diff_matrix));
        min_abs_diff = min_abs_diff(:);
        % - Identify correctly identified beats
        beat_correct_log = min_abs_diff<uParams.analysis.tol_window;
        clear diff_matrix min_abs_diff median_diff rel_ecg_els diff_matrix
        % - Calculate performance statistics
        sens(beat_detector_no, subj_no) = 100*sum(beat_correct_log)/length(beats_t(subj_no).ecg);
        ppv(beat_detector_no, subj_no) = 100*sum(beat_correct_log)/length(hq_ppg_beats);
        f1_score(beat_detector_no, subj_no) = 2*ppv(beat_detector_no, subj_no)*sens(beat_detector_no, subj_no)/(ppv(beat_detector_no, subj_no)+sens(beat_detector_no, subj_no));
        no_beats(beat_detector_no, subj_no) = length(beat_correct_log);
        clear hq_ppg_beats
        if do_acc
            % - Calculate performance according to accel
            min_val = max([beats_t(subj_no).ppg{beat_detector_no}(1)-rel_lag_t(subj_no), beats_t(subj_no).ecg(1)]);
            max_val = min([beats_t(subj_no).ppg{beat_detector_no}(end)-rel_lag_t(subj_no), beats_t(subj_no).ecg(end)]);
            win_durn = 20;
            win_starts = min_val:20:(max_val-win_durn);
            for win_no = 1 : length(win_starts)
                curr_win_start = win_starts(win_no);
                curr_win_end = win_starts(win_no)+win_durn;
                rel_ecg_beats = beats_t(subj_no).ecg>=curr_win_start & beats_t(subj_no).ecg <= curr_win_end;
                acc_counter = acc_counter+1;
                acc_wins.no_ecg{beat_detector_no}(acc_counter) = sum(rel_ecg_beats);
                acc_wins.no_ppg{beat_detector_no}(acc_counter) = sum((beats_t(subj_no).ppg{beat_detector_no}-rel_lag_t(subj_no)) >=curr_win_start & (beats_t(subj_no).ppg{beat_detector_no}-rel_lag_t(subj_no)) <=curr_win_end);
                acc_wins.no_correct{beat_detector_no}(acc_counter) = sum(beat_correct_log(rel_ecg_beats));
                rel_acc_els = ceil(curr_win_start*data(subj_no).acc_w.fs):(round(curr_win_end*data(subj_no).acc_w.fs)-1);
                acc_wins.mad{beat_detector_no}(acc_counter) = mean( abs( data(subj_no).acc_w.v(rel_acc_els) - mean(data(subj_no).acc_w.v(rel_acc_els)) ) );
            end
            clear win_no curr_win_start curr_win_end rel_ecg_beats rel_acc_els
        end
        clear beat_correct_log
        
        % --- Quality assessment
        
        
        % --- HR monitoring
        
        % - calculate intervals
        ppg_beats.t = beats_t(subj_no).ppg{beat_detector_no}-rel_lag_t(subj_no);
        ppg_beats.ibi = [diff(ppg_beats.t);nan];
        ecg_beats.t = beats_t(subj_no).ecg;
        ecg_beats.ibi = [diff(ecg_beats.t);nan];
        
        % - Generate time series
        ppg_int_ts = generate_beat_signal(ppg_beats.t, uParams);
        ecg_int_ts = generate_beat_signal(ecg_beats.t, uParams);
        
        % - Generate PPG quality time series
        ppg_beats.qual = beats_t(subj_no).ppg_qual{beat_detector_no};
        ppg_qual_ts = generate_qual_signal(ppg_beats.t, ppg_beats.qual, uParams);
    
        clear ppg_beats ecg_beats
        
        % - truncate to shortest time series
        min_t = max([ppg_qual_ts.t(find(~isnan(ppg_qual_ts.v), 1)), ppg_int_ts.t(find(~isnan(ppg_int_ts.v), 1)), ecg_int_ts.t(find(~isnan(ecg_int_ts.v), 1))]);
        max_t = min([ppg_qual_ts.t(find(~isnan(ppg_qual_ts.v), 1, 'last')), ppg_int_ts.t(find(~isnan(ppg_int_ts.v), 1, 'last')), ecg_int_ts.t(find(~isnan(ecg_int_ts.v), 1, 'last'))]);
        rel_els = ppg_int_ts.t >= min_t & ppg_int_ts.t<=max_t;
        ppg_int_ts.t = ppg_int_ts.t(rel_els);
        ppg_int_ts.v = ppg_int_ts.v(rel_els);
        rel_els = ppg_qual_ts.t >= min_t & ppg_qual_ts.t<=max_t;
        ppg_qual_ts.t = ppg_qual_ts.t(rel_els);
        ppg_qual_ts.v = ppg_qual_ts.v(rel_els);
        rel_els = ecg_int_ts.t >= min_t & ecg_int_ts.t<=max_t;
        ecg_int_ts.t = ecg_int_ts.t(rel_els);
        ecg_int_ts.v = ecg_int_ts.v(rel_els);
        clear min_t max_t rel_els
        
        % - eliminate times at which PPG qual was low
        total_durn = (ppg_int_ts.t(end)-ppg_int_ts.t(1)); % in secs
        rel_els = ppg_qual_ts.v==1;
        ppg_int_ts.t = ppg_int_ts.t(rel_els);
        ppg_int_ts.v = ppg_int_ts.v(rel_els);
        ecg_int_ts.t = ecg_int_ts.t(rel_els);
        ecg_int_ts.v = ecg_int_ts.v(rel_els);
        
        % - calculate the proportion of the time for which the PPG-dervied HR is within +/- 5bpm of the ECG-derived HR
        ppg_int_ts.hr = 60./ppg_int_ts.v;
        ecg_int_ts.hr = 60./ecg_int_ts.v;
        correct_hrs = ppg_int_ts.hr > (ecg_int_ts.hr-uParams.analysis.hr_tol) & ppg_int_ts.hr < (ecg_int_ts.hr+uParams.analysis.hr_tol);
        acc_hr(beat_detector_no, subj_no) = 100*mean(correct_hrs);
        durn_hr(beat_detector_no, subj_no) = (length(ppg_int_ts.t)/uParams.analysis.interpolation_fs)/60; % in mins
        prop_hr(beat_detector_no, subj_no) = 100*(length(ppg_int_ts.t)/uParams.analysis.interpolation_fs)/total_durn; % in percent
        bias_hr(beat_detector_no, subj_no) = mean(ppg_int_ts.hr - ecg_int_ts.hr);
        loa_hr(beat_detector_no, subj_no) = 1.96*std(ppg_int_ts.hr - ecg_int_ts.hr);
        
        clear correct_hrs ppg_int_ts ecg_int_ts
        
        
    end
    
    if do_acc
        % additional acc calcs
        acc_wins.sens{beat_detector_no} = 100*acc_wins.no_correct{beat_detector_no}./acc_wins.no_ecg{beat_detector_no};
        acc_wins.ppv{beat_detector_no} = 100*acc_wins.no_correct{beat_detector_no}./acc_wins.no_ppg{beat_detector_no};
        acc_wins.f1_score{beat_detector_no} = 2*acc_wins.ppv{beat_detector_no}.*acc_wins.sens{beat_detector_no}./(acc_wins.ppv{beat_detector_no}+acc_wins.sens{beat_detector_no});
    end
    
    % calculate results
    fprintf('\n   - Calculating results')
    stats = {'qual_sens', 'qual_spec', 'qual_bal_acc', 'sens', 'ppv', 'f1_score', 'no_beats', 'acc_hr', 'durn_hr', 'prop_hr', 'bias_hr', 'loa_hr'};
    summ_stats = {'med', 'lq', 'uq'};
    for stat_no = 1 : length(stats)
        curr_stat = stats{stat_no};
        eval(['curr_data = ' curr_stat '(beat_detector_no,:);']);
        curr_data = curr_data(~isnan(curr_data));
        
        for summ_stat_no = 1 : length(summ_stats)
            curr_summ_stat = summ_stats{summ_stat_no};
            switch curr_summ_stat
                case 'med'
                    temp = median(curr_data);
                case 'lq'
                    temp = quantile(curr_data, 0.25);
                case 'uq'
                    temp = quantile(curr_data, 0.75);
            end
            eval(['res.' curr_stat '.' curr_summ_stat '(beat_detector_no,1) = temp;']);
            clear temp curr_summ_stat
        end
        clear curr_stat curr_data summ_stat_no
    end
    clear stat_no
    
    % - output results
    fprintf('\n      -- Beat detector: %s', uParams.analysis.beat_detectors{beat_detector_no});
    for stat_no = 1 : length(stats)
        curr_stat = stats{stat_no};
        eval(['rel_res = res.' curr_stat ';']);
        fprintf('\n     -  %s: %.1f (%.1f - %.1f)', curr_stat, rel_res.med(beat_detector_no), rel_res.lq(beat_detector_no), rel_res.uq(beat_detector_no));
    end
    fprintf('\n     -  total beats: %d', sum(no_beats(beat_detector_no,:)));
    
    if do_acc
        % Create boxplot
        group_sep = 0.01;
        group_starts = 0:group_sep:0.3;
        box_data.labels = {};
        box_data.vals = [];
        for group_no = 1 : length(group_starts)
            curr_start = group_starts(group_no);
            curr_end = group_starts(group_no)+group_sep;
            rel_els = acc_wins.mad{beat_detector_no}(1:acc_counter) >= curr_start & ...
                acc_wins.mad{beat_detector_no}(1:acc_counter) < curr_end;
            no_correct = acc_wins.no_correct{beat_detector_no}(rel_els);
            no_ecg = acc_wins.no_ecg{beat_detector_no}(rel_els);
            no_ppg = acc_wins.no_ppg{beat_detector_no}(rel_els);
            ppv = sum(no_correct)/sum(no_ppg);
            sens = sum(no_correct)/sum(no_ecg);
            totals(group_no,1) = sum(no_ecg);
            all_vals(group_no,1) = 2*ppv*sens/(ppv+sens);
            box_data.vals(end+1:end+sum(rel_els)) = acc_wins.f1_score{beat_detector_no}(rel_els);
            box_data.labels(end+1:end+sum(rel_els)) = repmat({num2str(curr_start)}, sum(rel_els),1);
        end
        % boxplot(box_data.vals,box_data.labels)
        plot(all_vals)
    end
    
end


end

function assess_quality_assessment_tools(uParams)

fprintf('\n - Assessing performance of quality assessment tools: ')

%% See if there are any remaining quality tools to be assessed
% - create filepath
file_type = 'ppg_qual_perf';
filepath = create_proc_filepath(uParams, nan, file_type);
% - load file
if exist(filepath, 'file')
    load(filepath);  % loads 'ppg_qual_perf' variable
    quality_tools_used = fieldnames(ppg_qual_perf.raw);
    quality_tools_remaining = setxor(quality_tools_used, uParams.analysis.sig_qual_tools);
    % - skip this subject if there are no more beat detectors to use
    if isempty(quality_tools_remaining)
        fprintf('all done');
        return
    end
    clear ppg_qual_perf
end

% - go through all sig quality tools
quality_tools_remaining = uParams.analysis.sig_qual_tools;

% cycle through each subject
ppg_qual_perf.raw = struct;
for subj_no = 1 : uParams.dataset_details.no_subjs
    fprintf('\n - Assessing quality tools in subj %d: ', subj_no);
        
    %% Load time-aligned ECG beats
    file_type = 'ecg_beats_aligned';
    loadpath = create_proc_filepath(uParams, subj_no, file_type);
    load(loadpath); clear loadpath % loads ecg_beats_a_inds
    
    %% Load PPG beats
    file_type = 'ppg_beats';
    loadpath = create_proc_filepath(uParams, subj_no, file_type);
    load(loadpath); clear loadpath
    
    %% Load quality assessment results
    file_type = 'ppg_qual';
    loadpath = create_proc_filepath(uParams, subj_no, file_type);
    load(loadpath); clear loadpath
    
    %% Assess performance of each tool
    ppg_qual_perf.raw = assess_perf_quality_tools(ppg_qual, ppg_beats_inds, ecg_beats_a_inds, quality_tools_remaining, ppg_qual_perf.raw, subj_no, uParams);
    clear file_type ppg_beats_inds ecg_beats_a_inds ppg_qual
    
end
clear subj_no quality_tools_remaining 

%% Calculate summary statistics
quality_tools = fieldnames(ppg_qual_perf.raw);
for tool_no = 1 : length(quality_tools)
    curr_tool = quality_tools{tool_no};
    eval(['curr_tool_data = ppg_qual_perf.raw.' curr_tool ';']);
    perf_metrics = fieldnames(curr_tool_data);
    for metric_no = 1 : length(perf_metrics)
        curr_metric = perf_metrics{metric_no};
        eval(['ppg_qual_perf.stats.' curr_metric '.med(:,tool_no) = median(curr_tool_data.' curr_metric ',2);']);
        eval(['ppg_qual_perf.stats.' curr_metric '.lq(:,tool_no) = quantile(curr_tool_data.' curr_metric ',0.25,2);']);
        eval(['ppg_qual_perf.stats.' curr_metric '.uq(:,tool_no) = quantile(curr_tool_data.' curr_metric ',0.75,2);']);
    end
    clear curr_tool metric_no perf_metrics curr_tool_data curr_metric
end
ppg_qual_perf.stats.quality_tools = quality_tools;
ppg_qual_perf.stats.beat_detectors = uParams.analysis.beat_detectors;
clear tool_no quality_tools

%% Create ranking


%% Save quality tool performance results
save(filepath, 'ppg_qual_perf');

end

function ppg_qual_perf = assess_perf_quality_tools(ppg_qual, ppg_beats_inds, ecg_beats_a_inds, quality_tools_remaining, ppg_qual_perf, subj_no, uParams)

% - cycle through each quality assessment tool
for quality_tool_no = 1 : length(quality_tools_remaining)
    
    curr_tool = quality_tools_remaining{quality_tool_no};
    fprintf('%s, ', curr_tool);
    eval(['curr_ppg_qual = ppg_qual.' curr_tool '.v;']);
    
    for beat_detector_no = 1 : length(uParams.analysis.beat_detectors)
        
        curr_beat_detector = uParams.analysis.beat_detectors{beat_detector_no};
        
        % - create vector of ECG timings
        eval(['curr_ecg_beats_a_inds = ecg_beats_a_inds.' curr_beat_detector ';']);
        ecg_beats_t = (curr_ecg_beats_a_inds-1)./uParams.dataset_details.ecg_fs;
        
        % - PPG beat detections
        eval(['curr_ppg_beat_inds = ppg_beats_inds.' curr_beat_detector ';']);
        if isempty(curr_ppg_beat_inds)
            beat_correct_log = false(size(curr_ppg_beat_inds));
        else
            curr_ppg_beats_t = (curr_ppg_beat_inds-1)./uParams.dataset_details.ppg_fs_ds;
            % - difference matrix
            diff_matrix = repmat(curr_ppg_beats_t, [1, length(ecg_beats_t)]) - ecg_beats_t';
            % - Find minimum differences
            min_abs_diff = min(abs(diff_matrix'));
            min_abs_diff = min_abs_diff(:);
            % - Identify correct beats
            beat_correct_log = min_abs_diff<uParams.analysis.tol_window;
        end
        % - PPG quality values
        curr_ppg_beats_qual = false(size(curr_ppg_beat_inds));
        curr_ppg_beats_qual(curr_ppg_beat_inds>0) = curr_ppg_qual(curr_ppg_beat_inds(curr_ppg_beat_inds>0));
        % - performance statistics for this quality tool
        curr_ppg_qual_perf.spec(beat_detector_no, 1) = 100*sum(~beat_correct_log & ~curr_ppg_beats_qual)/sum(~beat_correct_log);
        curr_ppg_qual_perf.sens(beat_detector_no, 1) = 100*sum(beat_correct_log & curr_ppg_beats_qual)/sum(beat_correct_log);
        curr_ppg_qual_perf.bal_acc(beat_detector_no, 1) = 0.5*(curr_ppg_qual_perf.sens(beat_detector_no, 1)+curr_ppg_qual_perf.spec(beat_detector_no, 1));
        
        clear curr_ppg_beats_qual beat_correct_log min_abs_diff diff_matrix curr_ppg_beats_t curr_ppg_beat_inds curr_beat_detector
        
    end
    
    % store performance of this quality tool
    perf_metrics = {'sens', 'spec', 'bal_acc'};
    for metric_no = 1 : length(perf_metrics)
        curr_metric = perf_metrics{metric_no};
        eval(['ppg_qual_perf.' curr_tool '.' curr_metric '(:,subj_no) = curr_ppg_qual_perf.' curr_metric ';']);
    end
    
    clear beat_detector_no curr_ppg_qual_perf curr_tool
    
end

end

function assess_ppg_beat_detector_performance(uParams)

fprintf('\n - Assessing PPG beat detector performance')

%% See if there are any remaining beat detection strategies to be assessed
% - create filepath
file_type = 'ppg_detect_perf';
filepath = create_proc_filepath(uParams, nan, file_type);
all_strategies = create_beat_detection_strategies(uParams);
% - load file
if exist(filepath, 'file') && ~uParams.analysis.redo_analysis
    load(filepath);  % loads 'ppg_qual_perf' variable
    strategies_used = fieldnames(ppg_strategy_perf.raw);
    strategies_remaining = setxor(strategies_used, all_strategies);
    % - skip this subject if there are no more beat detectors to use
    if isempty(strategies_remaining)
        fprintf(': all done');
        return
    end
else
    strategies_remaining = all_strategies;
    ppg_strategy_perf.raw = struct;
end

% cycle through each subject
for subj_no = 1 : uParams.dataset_details.no_subjs
    fprintf('\n - Assessing PPG beat detection strategies on subj %d: ', subj_no);
        
    %% Load time-aligned PPG beats
    file_type = 'ppg_beats';
    loadpath = create_proc_filepath(uParams, subj_no, file_type);
    load(loadpath); clear loadpath % loads ppg_beats_inds and ppg_exc_log
    
    %% Load ECG beats
    file_type = 'ecg_beats_aligned';
    loadpath = create_proc_filepath(uParams, subj_no, file_type);
    load(loadpath); clear loadpath % loads ecg_beats_a_inds and ecg_exc_log
    
    %% Load quality assessment results
    file_type = 'ppg_qual';
    loadpath = create_proc_filepath(uParams, subj_no, file_type);
    load(loadpath); clear loadpath
    
    %% Assess performance of each tool
    ppg_strategy_perf.raw = assess_beat_detection_strategies(ppg_qual, ppg_beats_inds, ecg_beats_a_inds, ppg_exc_log, ecg_exc_a_log, strategies_remaining, ppg_strategy_perf.raw, subj_no, uParams);
    clear file_type ppg_beats_inds ecg_beats_a_inds ppg_qual
    
end
clear subj_no quality_tools_remaining 

%% Save performance results
save(filepath, 'ppg_strategy_perf');

end

function create_stats_and_tables(uParams)

fprintf('\n - Calculating stats and creating tables')

%% See if this has already been done
% - create filepath
file_type = 'ppg_detect_stats';
filepath = create_proc_filepath(uParams, nan, file_type);
all_strategies = create_beat_detection_strategies(uParams);
% - skip if this has been done
if exist(filepath, 'file') && ~uParams.analysis.redo_analysis
    fprintf(': all done');
    return
end

%% Load data
load_file_type = 'ppg_detect_perf';
load_filepath = create_proc_filepath(uParams, nan, load_file_type);
load(load_filepath, 'ppg_strategy_perf')

%% Calculate summary statistics
if ~sum(strcmp(fieldnames(ppg_strategy_perf), 'stats'))
    fprintf('\n    - Summary stats')
    all_strategies = fieldnames(ppg_strategy_perf.raw);
    for strategy_no = 1 : length(all_strategies)
        curr_tool = all_strategies{strategy_no};
        eval(['curr_strategy_data = ppg_strategy_perf.raw.' curr_tool ';']);
        perf_metrics = fieldnames(curr_strategy_data);
        for metric_no = 1 : length(perf_metrics)
            curr_metric = perf_metrics{metric_no};
            eval(['curr_data = curr_strategy_data.' curr_metric ';']);
            curr_data = curr_data(~isnan(curr_data));
            eval(['ppg_strategy_perf.stats.' curr_metric '.med(:,strategy_no) = median(curr_data,2);']);
            eval(['ppg_strategy_perf.stats.' curr_metric '.lq(:,strategy_no) = quantile(curr_data,0.25,2);']);
            eval(['ppg_strategy_perf.stats.' curr_metric '.uq(:,strategy_no) = quantile(curr_data,0.75,2);']);
            eval(['ppg_strategy_perf.stats.' curr_metric '.pc10(:,strategy_no) = quantile(curr_data,0.10,2);']);
            eval(['ppg_strategy_perf.stats.' curr_metric '.pc90(:,strategy_no) = quantile(curr_data,0.90,2);']);
            eval(['ppg_strategy_perf.stats.' curr_metric '.no_subjs(:,strategy_no) = length(curr_data);']);
        end
        clear curr_strategy metric_no perf_metrics curr_strategy_data curr_metric curr_data curr_tool
    end
    ppg_strategy_perf.stats.strategies = all_strategies;
    clear strategy_no
    
    % save
    save(load_filepath, 'ppg_strategy_perf');
end
clear load_filepath load_file_type

%% Create rankings
res = struct;

% - identify performance metrics
eval(['perf_metrics = fieldnames(ppg_strategy_perf.raw.' ppg_strategy_perf.stats.strategies{1} ');']);

% --- Best strategies with no quality assessment
fprintf('\n    - Rankings for strategies with no quality assessment')
rel_strategies = contains(ppg_strategy_perf.stats.strategies, '__none');
score_to_rank = 'f1_score';
ranking_direction = 'descend';
strategy_group_name = 'noQual';
res = create_table_of_top_ranked_strategies(ppg_strategy_perf, rel_strategies, strategy_group_name, score_to_rank, ranking_direction, res);
clear rel_strategies score_to_rank strategy_group_name

% --- Best strategies with accel only
fprintf('\n    - Rankings for strategies with accel only')
score_to_rank = 'acc_hr';
ranking_direction = 'descend';
rel_strategies = false(size(ppg_strategy_perf.stats.strategies));
for s = 1 : length(rel_strategies)
    if strcmp(ppg_strategy_perf.stats.strategies{s}(end-6:end), '__accel')
        rel_strategies(s) = true;
    end
end
if sum(rel_strategies) ~= 0      % - skip if there are no 'accel' strategies
    strategy_group_name = 'accel';
    res = create_table_of_top_ranked_strategies(ppg_strategy_perf, rel_strategies, strategy_group_name, score_to_rank, ranking_direction, res);
end
clear rel_strategies strategy_group_name

% --- Best strategies without accel with qual assessment using same beat detector
fprintf('\n    - Rankings for strategies with quality assessment and no accel')
%  - part 1 : identify best strategy for each beat detector
beat_dets = uParams.analysis.beat_detectors;
rel_strategies = false(size(ppg_strategy_perf.stats.strategies));
for beat_det_no = 1 : length(beat_dets)
    curr_beat_det = beat_dets{beat_det_no};
    temp_rel_strategies = false(size(ppg_strategy_perf.stats.strategies));
    for s = 1 : length(temp_rel_strategies)
        if ~contains(ppg_strategy_perf.stats.strategies{s}, 'accel') && length(strfind(ppg_strategy_perf.stats.strategies{s}, curr_beat_det))==2
            temp_rel_strategies(s) = true;
        end
    end
    temp_rel_strategies = find(temp_rel_strategies);
    if ~isempty(temp_rel_strategies)
        eval(['curr_perf_stats = ppg_strategy_perf.stats.' score_to_rank '.med(temp_rel_strategies);']);
        [~,order] = sort(curr_perf_stats, 'descend');
        rel_strategies(temp_rel_strategies(order(1))) = true;
    end
end
%  - part 2 : create table of these strategies
if sum(rel_strategies) ~= 0      % - skip if there are no 'accel' strategies
    strategy_group_name = 'comb_noAccel';
    res = create_table_of_top_ranked_strategies(ppg_strategy_perf, rel_strategies, strategy_group_name, score_to_rank, ranking_direction, res);
end
clear rel_strategies strategy_group_name


% --- Best strategies with accel with qual assessment using same beat detector
fprintf('\n    - Rankings for strategies with quality assessment and accel')
%  - part 1 : identify best strategy for each beat detector
beat_dets = uParams.analysis.beat_detectors;
rel_strategies = false(size(ppg_strategy_perf.stats.strategies));
for beat_det_no = 1 : length(beat_dets)
    curr_beat_det = beat_dets{beat_det_no};
    temp_rel_strategies = false(size(ppg_strategy_perf.stats.strategies));
    for s = 1 : length(temp_rel_strategies)
        if contains(ppg_strategy_perf.stats.strategies{s}, 'accel') && length(strfind(ppg_strategy_perf.stats.strategies{s}, curr_beat_det))==2
            temp_rel_strategies(s) = true;
        end
    end
    temp_rel_strategies = find(temp_rel_strategies);
    eval(['curr_perf_stats = ppg_strategy_perf.stats.' score_to_rank '.med(temp_rel_strategies);']);
    [~,order] = sort(curr_perf_stats, 'descend');
    if ~isempty(order)
        rel_strategies(temp_rel_strategies(order(1))) = true;
    end
end
%  - part 2 : create table of these strategies
if sum(rel_strategies) ~= 0      % - skip if there are no 'accel' strategies
    strategy_group_name = 'comb_Accel';
    res = create_table_of_top_ranked_strategies(ppg_strategy_perf, rel_strategies, strategy_group_name, score_to_rank, ranking_direction, res);
end


%% Save performance results
save(filepath, 'res');

%% Remaining tables to be implemented

return


[~,order] = sort(ppg_strategy_perf.stats.new_score.med(rel_strategies), 'descend');
els = find(rel_strategies); els = els(order);
for metric_no = 1 : length(perf_metrics)
    curr_metric = perf_metrics{metric_no};
    for el_no = 1 : length(els)
        eval(['med = ppg_strategy_perf.stats.' curr_metric '.med(els(el_no));']);
        eval(['lq = ppg_strategy_perf.stats.' curr_metric '.lq(els(el_no));']);
        eval(['uq = ppg_strategy_perf.stats.' curr_metric '.uq(els(el_no));']);
        curr_txt = sprintf('%.1f (%.1f - %.1f)', med, lq, uq);
        eval([curr_metric '{el_no,1} = curr_txt;']);
        % -  name of this strategy (i.e. beat detector)
        eval(['strategy{el_no,1} = ppg_strategy_perf.stats.strategies(els(el_no));']);
        strategy{el_no} = strategy{el_no}; clear temp
        clear med lq uq curr_txt
    end
    clear el_no curr_metric
end
clear metric_no els
res.comb_no_accel = table(strategy, sens, ppv, f1_score, new_score);
clear strategy sens ppv f1_score new_score


% --- Best strategies without accel with qual assessment using any beat detector
%  - part 1 : identify best strategy for each beat detector
beat_dets = uParams.analysis.beat_detectors;
rel_strategies = false(size(ppg_strategy_perf.stats.strategies));
for beat_det_no = 1 : length(beat_dets)
    curr_beat_det = beat_dets{beat_det_no};
    temp_rel_strategies = false(size(ppg_strategy_perf.stats.strategies));
    for s = 1 : length(temp_rel_strategies)
        temp = strfind(ppg_strategy_perf.stats.strategies{s}, '__');
        beat_detector_name = ppg_strategy_perf.stats.strategies{s}(1:temp(1)-1);
        clear temp
        if ~contains(ppg_strategy_perf.stats.strategies{s}, 'accel') && strcmp(beat_detector_name, curr_beat_det)
            temp_rel_strategies(s) = true;
        end
    end
    temp_rel_strategies = find(temp_rel_strategies);
    [~,order] = sort(ppg_strategy_perf.stats.ppv.med(temp_rel_strategies), 'descend');
    rel_strategies(temp_rel_strategies(order(1))) = true;
end
%  - part 2 : sort these strategies
[~,order] = sort(ppg_strategy_perf.stats.new_score.med(rel_strategies), 'descend');
els = find(rel_strategies); els = els(order);
for metric_no = 1 : length(perf_metrics)
    curr_metric = perf_metrics{metric_no};
    for el_no = 1 : length(els)
        eval(['med = ppg_strategy_perf.stats.' curr_metric '.med(els(el_no));']);
        eval(['lq = ppg_strategy_perf.stats.' curr_metric '.lq(els(el_no));']);
        eval(['uq = ppg_strategy_perf.stats.' curr_metric '.uq(els(el_no));']);
        curr_txt = sprintf('%.1f (%.1f - %.1f)', med, lq, uq);
        eval([curr_metric '{el_no,1} = curr_txt;']);
        % -  name of this strategy (i.e. beat detector)
        eval(['strategy{el_no,1} = ppg_strategy_perf.stats.strategies(els(el_no));']);
        clear med lq uq curr_txt
    end
    clear el_no curr_metric
end
clear metric_no els
res.comb_any_no_accel = table(strategy, sens, ppv, f1_score, new_score);
clear strategy sens ppv f1_score new_score

end

function res = create_table_of_top_ranked_strategies(ppg_strategy_perf, rel_strategies, strategy_group_name, score_to_rank, ranking_direction, res)
eval(['scores_for_ranking = ppg_strategy_perf.stats.' score_to_rank '.med;']);
[~,order] = sort(scores_for_ranking(rel_strategies), ranking_direction);
els = find(rel_strategies); els = els(order);
eval(['perf_metrics = fieldnames(ppg_strategy_perf.raw.' ppg_strategy_perf.stats.strategies{find(rel_strategies,1)} ');']);
stat_types = {'med', 'pc10', 'lq', 'uq', 'pc90'};
for metric_no = 1 : length(perf_metrics)
    curr_metric = perf_metrics{metric_no};
    for el_no = 1 : length(els)
        for stat_type_no = 1:length(stat_types)
            curr_stat = stat_types{stat_type_no};
            eval([curr_stat ' = ppg_strategy_perf.stats.' curr_metric '.' curr_stat '(els(el_no));']);
        end
        curr_txt = sprintf('%.1f (%.1f - %.1f)', med, lq, uq);
        eval(['metric_res.' curr_metric '.txt{el_no,1} = curr_txt;']);
        for stat_type_no = 1:length(stat_types)
            curr_stat = stat_types{stat_type_no};
            eval(['metric_res.' curr_metric '.num.' curr_stat '(el_no,1) = ' curr_stat ';']);
        end
        % -  name of this strategy (i.e. beat detector)
        strategy_text = ppg_strategy_perf.stats.strategies{els(el_no)};
        temp = strfind(strategy_text, '__');
        strategy{el_no,1} = strategy_text(1:temp(1)-1);
        % - name of second beat detector (used for quality assessment)
        if strcmp(strategy_text(temp(1):end), '__none')
            qual{el_no,1} = 'none';
        elseif strcmp(strategy_text(temp(1):end), '__accel')
            qual{el_no,1} = 'accel';
        else
            temp = strfind(strategy_text, 'comb_');
            temp2 = strfind(strategy_text, '_'); temp2 = temp2(temp2>(temp(1)));
            qual_beat_detectors{1} = strategy_text(temp2(1)+1:temp2(2)-1);
            qual_beat_detectors{2} = strategy_text(temp2(2)+1:end);
            qual{el_no,1} = qual_beat_detectors{~strcmp(qual_beat_detectors,strategy{el_no})};
        end
        clear temp strategy_text
        clear med lq uq curr_txt
    end
    clear el_no curr_metric
end
clear metric_no els
% Create text table
cols = [];
for metric_no = 1 : length(perf_metrics)
    eval([perf_metrics{metric_no} ' = metric_res.' perf_metrics{metric_no} '.txt;']);
    cols = [cols, perf_metrics{metric_no}, ', '];
end
cols = cols(1:end-2);
try
    eval(['res.' strategy_group_name '.txt = table(strategy, qual, ' cols ');']);
catch
    a = 1;
end
clear cols metric_no
% Create numerical table
cols = [];
for metric_no = 1 : length(perf_metrics)
    for stat_type_no = 1:length(stat_types)
        curr_stat = stat_types{stat_type_no};
        eval([perf_metrics{metric_no} '_' curr_stat ' = metric_res.' perf_metrics{metric_no} '.num.' curr_stat ';']);
        cols = [cols, perf_metrics{metric_no} '_' curr_stat , ', '];
    end
end
cols = cols(1:end-2);
eval(['res.' strategy_group_name '.num = table(strategy, ' cols ');']);
clear cols metric_no
clear strategy sens ppv f1_score
end

function ppg_strategy_perf = assess_beat_detection_strategies(ppg_qual, ppg_beats_inds, ecg_beats_a_inds, ppg_exc_log, ecg_exc_a_log, strategies_remaining, ppg_strategy_perf, subj_no, uParams)

% - insert nans if there are no suitable data for analysis:
if unique(ppg_exc_log) == 1
    perf_metrics = {'sens'; 'ppv'; 'f1_score'; 'no_beats'; 'acc_hr'; 'durn_hr'; 'prop_hr'; 'bias_hr'; 'loa_hr'; 'mae_hr'; 'acc_ibi'; 'durn_ibi'; 'prop_ibi'; 'bias_ibi'; 'loa_ibi'; 'mae_ibi'};
    for strategy_no = 1 : length(strategies_remaining)
        curr_strategy = strategies_remaining{strategy_no};
        for metric_no = 1 : length(perf_metrics)
            curr_metric = perf_metrics{metric_no};
            eval(['ppg_strategy_perf.' curr_strategy '.' curr_metric '(1,subj_no) = nan;']);
        end
    end
    return
end

% - cycle through each beat detection strategy
for strategy_no = 1 : length(strategies_remaining)
    
    % - extract details of this strategy
    curr_strategy = strategies_remaining{strategy_no};
    temp = strfind(curr_strategy, '__');
    curr_beat_detector = curr_strategy(1:temp(1)-1);
    curr_qual_tool = curr_strategy(temp(1)+2:end);
    clear temp
    
    % --- extract required data ---
    
    % - extract data relating to this strategy
    eval(['curr_ppg_beat_inds = ppg_beats_inds.' curr_beat_detector ';']);
    eval(['curr_ecg_beat_a_inds = ecg_beats_a_inds.' curr_beat_detector ';']);
    eval(['curr_ppg_qual = ppg_qual.' curr_qual_tool ';']);
    clear curr_qual_tool curr_beat_detector
    
    % - evaluate quality of each PPG beat detection
    ppg_beat_qual_log = create_ppg_beat_qual_log(curr_ppg_beat_inds, curr_ppg_qual);
    
    % - create vectors of PPG and ECG beat detection timings
    curr_ppg_beats.t = (curr_ppg_beat_inds-1)./uParams.dataset_details.ppg_fs_ds;
    curr_ecg_beats.t = (curr_ecg_beat_a_inds-1)./uParams.dataset_details.ecg_fs;
    clear curr_ecg_beat_a_inds curr_ppg_beat_inds curr_ppg_qual
    
    % - calculate inter-beat intervals
    curr_ppg_beats.ibi = [diff(curr_ppg_beats.t);nan];
    curr_ecg_beats.ibi = [diff(curr_ecg_beats.t);nan];
    
    % - create logicals indicating whether or not beats should be included in the assessment
    curr_ppg_beats.inc_log = create_inclusion_logical(curr_ppg_beats.t, ppg_exc_log, uParams.dataset_details.ppg_fs_ds, ecg_exc_a_log, uParams.dataset_details.ecg_fs);
    curr_ecg_beats.inc_log = create_inclusion_logical(curr_ecg_beats.t, ppg_exc_log, uParams.dataset_details.ppg_fs_ds, ecg_exc_a_log, uParams.dataset_details.ecg_fs);
    
    % - refine the vectors of PPG and ECG beat detection timings to only include those which are to be included in the analysis
    curr_ppg_beats_inc.t = curr_ppg_beats.t(curr_ppg_beats.inc_log);
    curr_ppg_beats_inc.ibi = curr_ppg_beats.ibi(curr_ppg_beats.inc_log);
    ppg_beat_qual_log_inc = ppg_beat_qual_log(curr_ppg_beats.inc_log);
    curr_ecg_beats_inc.t = curr_ecg_beats.t(curr_ecg_beats.inc_log);
    curr_ecg_beats_inc.ibi = curr_ecg_beats.ibi(curr_ecg_beats.inc_log);
    
    % --- beat detection ---
    
    % - Calculate time differences between each pair of PPG and ECG beats
    hq_ppg_beats = curr_ppg_beats_inc.t(ppg_beat_qual_log_inc);
    if isempty(hq_ppg_beats)
        beat_correct_log = false(size(curr_ecg_beats_inc.t));
    else
        diff_matrix = repmat(hq_ppg_beats, [1, length(curr_ecg_beats_inc.t)]) - curr_ecg_beats_inc.t';
        % - Find minimum differences
        min_abs_diff = min(abs(diff_matrix), [], 2); % used to be: min_abs_diff = min(abs(diff_matrix));
        %min_abs_diff = min_abs_diff(:);
        % - Identify correctly identified beats
        beat_correct_log = min_abs_diff<uParams.analysis.tol_window;
    end
    clear diff_matrix min_abs_diff
    % - Calculate performance statistics
    curr_ppg_strategy_perf.sens = 100*sum(beat_correct_log)/length(curr_ecg_beats_inc.t);
    curr_ppg_strategy_perf.ppv = 100*sum(beat_correct_log)/length(hq_ppg_beats);
    curr_ppg_strategy_perf.f1_score = 2*curr_ppg_strategy_perf.ppv*curr_ppg_strategy_perf.sens/(curr_ppg_strategy_perf.ppv+curr_ppg_strategy_perf.sens);
    %curr_ppg_strategy_perf.new_score = (5*curr_ppg_strategy_perf.ppv+curr_ppg_strategy_perf.sens)/6;
    curr_ppg_strategy_perf.no_beats = length(curr_ecg_beats_inc.t);
    clear hq_ppg_beats
    
    clear beat_detector_no curr_tool beat_correct_log hq_ppg_beats ppg_beat_qual_log_inc curr_ecg_beats_inc curr_ppg_beats_inc
    
    % --- HR estimation
    
    % - Calculate HRs using windowing
    signal_start_t = 0;
    signal_end_t = max(curr_ppg_beats.t);
    curr_ecg_beats.hr = calc_hrs_from_beat_timings(curr_ecg_beats.t, curr_ecg_beats.inc_log, signal_start_t, signal_end_t, uParams.analysis.hr_win_durn);
    beat_qual = ppg_beat_qual_log & curr_ppg_beats.inc_log;
    curr_ppg_beats.hr = calc_hrs_from_beat_timings(curr_ppg_beats.t, beat_qual, signal_start_t, signal_end_t, uParams.analysis.hr_win_durn);
    clear beat_qual
    % CHECK:   plot(curr_ecg_beats.t, curr_ecg_beats.hr), hold on, plot(curr_ppg_beats.t, curr_ppg_beats.hr)
    
    % - Generate time series of HRs
    ppg_hr_ts = generate_hr_signal(curr_ppg_beats.t, curr_ppg_beats.hr, signal_start_t, signal_end_t, uParams);
    ecg_hr_ts = generate_hr_signal(curr_ecg_beats.t, curr_ecg_beats.hr, signal_start_t, signal_end_t, uParams);
    clear signal_start_t signal_end_t
    
    % - Generate PPG quality time series
    ppg_qual_and_inc_ts = generate_qual_signal(curr_ppg_beats.t, ppg_beat_qual_log & curr_ppg_beats.inc_log, uParams);
    
    % - Generate moving window PPG quality time series
    no_samps = ceil(uParams.analysis.hr_win_durn*uParams.analysis.interpolation_fs);
    ppg_qual_and_inc_ts.v_win = movmin(ppg_qual_and_inc_ts.v,no_samps);
    
    % - truncate to shortest time series
    min_t = max([ppg_qual_and_inc_ts.t(find(~isnan(ppg_qual_and_inc_ts.v_win), 1)), ppg_hr_ts.t(find(~isnan(ppg_hr_ts.v), 1)), ecg_hr_ts.t(find(~isnan(ecg_hr_ts.v), 1))]);
    max_t = min([ppg_qual_and_inc_ts.t(find(~isnan(ppg_qual_and_inc_ts.v_win), 1, 'last')), ppg_hr_ts.t(find(~isnan(ppg_hr_ts.v), 1, 'last')), ecg_hr_ts.t(find(~isnan(ecg_hr_ts.v), 1, 'last'))]);
    rel_els = ppg_hr_ts.t >= min_t & ppg_hr_ts.t<=max_t;
    ppg_hr_ts.t = ppg_hr_ts.t(rel_els);
    ppg_hr_ts.v = ppg_hr_ts.v(rel_els);
    rel_els = ppg_qual_and_inc_ts.t >= min_t & ppg_qual_and_inc_ts.t<=max_t;
    ppg_qual_and_inc_ts.t = ppg_qual_and_inc_ts.t(rel_els);
    ppg_qual_and_inc_ts.v = ppg_qual_and_inc_ts.v(rel_els);
    ppg_qual_and_inc_ts.v_win = ppg_qual_and_inc_ts.v_win(rel_els);
    rel_els = ecg_hr_ts.t >= min_t & ecg_hr_ts.t<=max_t;
    ecg_hr_ts.t = ecg_hr_ts.t(rel_els);
    ecg_hr_ts.v = ecg_hr_ts.v(rel_els);
    clear min_t max_t rel_els
    
    if isempty(ppg_hr_ts.t)
        % - fill in results if there are no data available
        [a.acc_hr, a.durn_hr, a.prop_hr, a.bias_hr, a.loa_hr, a.mae_hr] = deal(nan);
        curr_ppg_strategy_perf = a; clear a
    else
        % - eliminate times at which PPG qual was low
        total_durn = (ppg_hr_ts.t(end)-ppg_hr_ts.t(1))/60; % in mins
        rel_els = ppg_qual_and_inc_ts.v_win==1;
        ppg_hr_ts.t = ppg_hr_ts.t(rel_els);
        ppg_hr_ts.v = ppg_hr_ts.v(rel_els);
        ecg_hr_ts.t = ecg_hr_ts.t(rel_els);
        ecg_hr_ts.v = ecg_hr_ts.v(rel_els);
        clear rel_els ppg_qual_and_inc_ts
        
        % - calculate the proportion of the time for which the PPG-dervied HR is within +/- 5bpm of the ECG-derived HR
        hr_errors = ppg_hr_ts.v - ecg_hr_ts.v;
        correct_hrs = ppg_hr_ts.v > (ecg_hr_ts.v-uParams.analysis.hr_tol) & ppg_hr_ts.v < (ecg_hr_ts.v+uParams.analysis.hr_tol);
        curr_ppg_strategy_perf.acc_hr = 100*mean(correct_hrs);
        curr_ppg_strategy_perf.durn_hr = (length(ppg_hr_ts.t)/uParams.analysis.interpolation_fs)/60; % in mins
        curr_ppg_strategy_perf.prop_hr = 100*(curr_ppg_strategy_perf.durn_hr)/total_durn; % in percent
        curr_ppg_strategy_perf.bias_hr = mean(hr_errors);
        curr_ppg_strategy_perf.loa_hr = 1.96*std(hr_errors);
        curr_ppg_strategy_perf.mae_hr = mean(abs(hr_errors));
        
    end
    
    clear correct_hrs ppg_hr_ts ecg_hr_ts total_durn hr_errors
    
    % store performance of this beat detection strategy
    perf_metrics = fieldnames(curr_ppg_strategy_perf);
    for metric_no = 1 : length(perf_metrics)
        curr_metric = perf_metrics{metric_no};
        eval(['ppg_strategy_perf.' curr_strategy '.' curr_metric '(1,subj_no) = curr_ppg_strategy_perf.' curr_metric ';']);
    end
    clear metric_no perf_metrics curr_metric curr_ppg_strategy_perf
    
    
    % --- IBI measurement
    
    % - Generate time series
    ppg_int_ts = generate_beat_signal(curr_ppg_beats.t, uParams);
    ecg_int_ts = generate_beat_signal(curr_ecg_beats.t, uParams);
    
    % - Generate PPG quality time series
    ppg_qual_and_inc_ts = generate_qual_signal(curr_ppg_beats.t, ppg_beat_qual_log & curr_ppg_beats.inc_log, uParams);
    
    clear ppg_beats ecg_beats ppg_beat_qual_log curr_ecg_beats curr_ppg_beats
    
    % - truncate to shortest time series
    min_t = max([ppg_qual_and_inc_ts.t(find(~isnan(ppg_qual_and_inc_ts.v), 1)), ppg_int_ts.t(find(~isnan(ppg_int_ts.v), 1)), ecg_int_ts.t(find(~isnan(ecg_int_ts.v), 1))]);
    max_t = min([ppg_qual_and_inc_ts.t(find(~isnan(ppg_qual_and_inc_ts.v), 1, 'last')), ppg_int_ts.t(find(~isnan(ppg_int_ts.v), 1, 'last')), ecg_int_ts.t(find(~isnan(ecg_int_ts.v), 1, 'last'))]);
    rel_els = ppg_int_ts.t >= min_t & ppg_int_ts.t<=max_t;
    ppg_int_ts.t = ppg_int_ts.t(rel_els);
    ppg_int_ts.v = ppg_int_ts.v(rel_els);
    rel_els = ppg_qual_and_inc_ts.t >= min_t & ppg_qual_and_inc_ts.t<=max_t;
    ppg_qual_and_inc_ts.t = ppg_qual_and_inc_ts.t(rel_els);
    ppg_qual_and_inc_ts.v = ppg_qual_and_inc_ts.v(rel_els);
    rel_els = ecg_int_ts.t >= min_t & ecg_int_ts.t<=max_t;
    ecg_int_ts.t = ecg_int_ts.t(rel_els);
    ecg_int_ts.v = ecg_int_ts.v(rel_els);
    clear min_t max_t rel_els
    
    % - eliminate times at which PPG qual was low
    total_durn = (ppg_int_ts.t(end)-ppg_int_ts.t(1))/60; % in mins
    rel_els = ppg_qual_and_inc_ts.v==1;
    ppg_int_ts.t = ppg_int_ts.t(rel_els);
    ppg_int_ts.v = ppg_int_ts.v(rel_els);
    ecg_int_ts.t = ecg_int_ts.t(rel_els);
    ecg_int_ts.v = ecg_int_ts.v(rel_els);
    clear rel_els ppg_qual_and_inc_ts
    
    % - calculate the proportion of the time for which the PPG-dervied HR is within +/- 5bpm of the ECG-derived HR
    ppg_int_ts.hr = 60./ppg_int_ts.v;
    ecg_int_ts.hr = 60./ecg_int_ts.v;
    ibi_errors = 1000*(ppg_int_ts.v - ecg_int_ts.v);
    correct_ibis = ppg_int_ts.hr > (ecg_int_ts.hr-uParams.analysis.hr_tol) & ppg_int_ts.hr < (ecg_int_ts.hr+uParams.analysis.hr_tol);
    curr_ppg_strategy_perf.acc_ibi = 100*mean(correct_ibis);
    curr_ppg_strategy_perf.durn_ibi = (length(ppg_int_ts.t)/uParams.analysis.interpolation_fs)/60; % in mins
    curr_ppg_strategy_perf.prop_ibi = 100*(curr_ppg_strategy_perf.durn_ibi)/total_durn; % in percent
    curr_ppg_strategy_perf.bias_ibi = mean(ibi_errors);
    curr_ppg_strategy_perf.loa_ibi = 1.96*std(ibi_errors);
    curr_ppg_strategy_perf.mae_ibi = mean(abs(ibi_errors));
    
    clear correct_ibis ppg_int_ts ecg_int_ts total_durn 
    
    % store performance of this beat detection strategy
    perf_metrics = fieldnames(curr_ppg_strategy_perf);
    for metric_no = 1 : length(perf_metrics)
        curr_metric = perf_metrics{metric_no};
        eval(['ppg_strategy_perf.' curr_strategy '.' curr_metric '(1,subj_no) = curr_ppg_strategy_perf.' curr_metric ';']);
    end
    clear metric_no perf_metrics curr_strategy curr_metric curr_ppg_strategy_perf
    
end

end

function hr = calc_hrs_from_beat_timings(beats_t, beat_qual, signal_start_t, signal_end_t, win_durn)

hr = nan(size(beats_t));
for beat_no = 1 : length(hr)
    % - calculate timings of this window for HR estimation
    end_t = beats_t(beat_no);
    start_t = end_t - win_durn;
    % skip if there isn't complete data for this window
    if start_t < signal_start_t || end_t > signal_end_t
        continue
    end
    % - identify beats in this window
    rel_beats = beats_t >= start_t & beats_t <= end_t;
    % skip if there is low quality PPG data in this window
    if sum(~beat_qual(rel_beats))
        continue
    end
    % - calculate HRs
    hr(beat_no) = calc_hr_from_window(beats_t(rel_beats), win_durn);
    
end



end

function hr = calc_hr_from_window(beats_t, win_durn)

% skip if no beats were detected in this window
if isempty(beats_t)
    hr = nan;
    return
end

% calculate HR
no_beats = length(beats_t)-1; % take off one to give the number of complete heart beats in this time period
time_diff = range(beats_t); % in secs
hr = 60*(no_beats/time_diff); % in bpm

end

function inc_log = create_inclusion_logical(curr_beats_t, ppg_exc_log, ppg_fs, ecg_exc_log, ecg_fs)

min_t = 0;
max_t = min([(length(ppg_exc_log)-1)./ppg_fs, (length(ecg_exc_log)-1)./ecg_fs]);
rel_inds = curr_beats_t>min_t & curr_beats_t<max_t;

curr_beats_ppg_inds = round(curr_beats_t(rel_inds)*ppg_fs) + 1;
rel_ppg_exc_log = true(length(curr_beats_t),1);
rel_ppg_exc_log(rel_inds) = ppg_exc_log(curr_beats_ppg_inds);

curr_beats_ecg_inds = round(curr_beats_t(rel_inds)*ecg_fs) + 1;
rel_ecg_exc_log = true(length(curr_beats_t),1);
rel_ecg_exc_log(rel_inds) = ecg_exc_log(curr_beats_ecg_inds);

inc_log = ~rel_ppg_exc_log & ~rel_ecg_exc_log;

end

function ppg_beat_qual_log = create_ppg_beat_qual_log(curr_ppg_beat_inds, curr_ppg_qual)
% create a logical indicating whether or not each ppg beat detection was during a period of high or low quality.

ppg_beat_qual_log = curr_ppg_qual.v(curr_ppg_beat_inds);

end

function strategies = create_beat_detection_strategies(uParams)

% - specify beat detectors
beat_detectors = uParams.analysis.beat_detectors;

% - specify signal quality tools
qual_tools = uParams.analysis.sig_qual_tools;

% - create list of beat detection strategies
strategies = cell(length(beat_detectors)*length(qual_tools),1);
counter_no = 0;
for det_no = 1 : length(beat_detectors)
    curr_det = beat_detectors{det_no};
    
    for tool_no = 1 : length(qual_tools)
        curr_tool = qual_tools{tool_no};
        
        counter_no = counter_no+1;
        strategies{counter_no} = [curr_det, '__', curr_tool];
    end
    
end

end

function uParams = setup_up(dataset, options)
% Sets up universal parameters for use throughout the code

%% Display startup message
display_startup_message;

%% File path of dataset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% FILE PATH TO BE SPECIFIED %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% (may require editing) %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% - If the path of the dataset file has been specified in the 'options' structure, then copy it across:
if sum(strcmp(fieldnames(options), 'dataset_file'))
    uParams.paths.dataset_file = options.dataset_file;
else
    % Otherwise, the following hard-coded path is used for the dataset file:
    dataset_filename = 'unknown_data.mat';          % This should contain the name of the dataset file
    root_folder = 'C:/unknown_dataset_folder/'; % This should contain the path of the directory containing the dataset file
    uParams.paths.dataset_file = [root_folder, dataset_filename];
    
    fprintf('\n - Using the following hard-coded path:\n   %s', uParams.paths.dataset_file);
    options.dataset_file = dataset_filename; clear dataset_filename
end

%% Name of dataset
if isempty(dataset)
    uParams.dataset = strrep(options.dataset_file, '_data', '');
    uParams.dataset = strrep(uParams.dataset, '.mat', '');
else
    uParams.dataset = dataset; clear dataset
end

%% Display name of dataset being processed
fprintf('\n ~~~ PPG Beat Detector Assessment ~~~')
fprintf(['\n  ~~~ on ', uParams.dataset, ' Dataset ~~~\n'])

%% Setup default options
uParams = setup_default_options(uParams,options);

%% Other paths
[folder_containing_raw_data_file, ~, ~] = fileparts(options.dataset_file);
eval(sprintf('uParams.paths.processing_folder = ''%s%s%s%s%s'';', folder_containing_raw_data_file, filesep, 'proc_data_', uParams.dataset, filesep));
if ~exist(uParams.paths.processing_folder, 'dir')
    fprintf('\n - Making directory to store results in at: %s', uParams.paths.processing_folder);
    mkdir(uParams.paths.processing_folder);
end

%% Analysis parameters
uParams.analysis.tol_window = 0.2; % in secs
uParams.analysis.qrs_tol_window = 0.2; % in secs
uParams.analysis.win_durn = 20; % in secs
uParams.analysis.win_overlap = 5; % in secs
uParams.analysis.ppg_fid_pt = 'mid_pt'; % options: mid_pt, pk, onset
uParams.analysis.do_wins = 1; % whether or not to analyse a PPG signal in windows rather than as one long segment
uParams.analysis.interpolation_fs = 50; % sampling freq of interpolated inter-beat-interval time series (shouldn't be more than 100 Hz)
uParams.analysis.max_lag = 10; % max permissible lag between ECG and PPG in secs
uParams.analysis.lag_int = 0.02;
uParams.analysis.hr_tol = 5; % tolerance in bpm of HR estimates (to be classified as accurate)
uParams.analysis.durn_flat_line = 0.1; % threshold duration in secs, above which a flat line is considered to indicate no signal (below this it might just be due to a temporarily steady PPG value).
uParams.analysis.sig_qual_tools = {}; % {'accel'; 'comb_beat_detectors'};  % a list of signal quality assessment tools to use
uParams.analysis.sig_qual_tools = add_in_comb_beat_detectors(uParams.analysis.sig_qual_tools, uParams);
uParams.analysis.hr_win_durn = 8; % Duration of window (in secs) over which to estimate HRs.

%% Combined PPG beat detector quality assessment
uParams.analysis.window_durn_for_incorrect_beats = 5; % in secs

%% Filtering parameters

% - Fiducial Points: eliminating high freqs
uParams.analysis.filtering.elim_high_freqs.Fpass = 10;   % in HZ
uParams.analysis.filtering.elim_high_freqs.Fstop = 6.93;   % in HZ   % 6.93 - 10 gives -3dB cutoff of 8.00 Hz
uParams.analysis.filtering.elim_high_freqs.Dpass = 0.05;
uParams.analysis.filtering.elim_high_freqs.Dstop = 0.01;

% - Fiducial Points: eliminating low freqs
uParams.analysis.filtering.elim_low_freqs.Fpass = 1.02;  % in Hz
uParams.analysis.filtering.elim_low_freqs.Fstop = 0.6;  % in Hz      % 0.6 - 1.02 gives -3dB cutoff of 40 bpm (i.e. 0.67 Hz)
uParams.analysis.filtering.elim_low_freqs.Dpass = 0.05;
uParams.analysis.filtering.elim_low_freqs.Dstop = 0.01;

%% Quality assessment parameters
uParams.qual_assessment.do_acc = 1;  % whether or not to assess quality using a simultaneous accelerometry signal

%% Dataset-specific settings
% - Find out how many subjects are in the dataset
data = load_data(uParams);
uParams.dataset_details.no_subjs = length(data);
uParams.dataset_details.ppg_fs = data(1).ppg.fs;
if uParams.analysis.do_downsample && (uParams.analysis.downsample_freq < data(1).ppg.fs)
    uParams.dataset_details.ppg_fs_ds = uParams.analysis.downsample_freq;
else
    uParams.analysis.do_downsample = 0;
    uParams.dataset_details.ppg_fs_ds = uParams.dataset_details.ppg_fs;
end
uParams.dataset_details.ecg_fs = data(1).ecg.fs;
% - Remove any signal quality assessment tools which can't be used on this dataset
%   - accel
if ~sum(strcmp(fieldnames(data), 'acc_ppg_site'))
    uParams.analysis.sig_qual_tools = uParams.analysis.sig_qual_tools(~contains(uParams.analysis.sig_qual_tools, 'accel'));
end

end

function sig_qual_tools = add_in_comb_beat_detectors(sig_qual_tools, uParams)

% - add in 'none' - i.e. not using quality assessment
sig_qual_tools = ['none'; sig_qual_tools(:)];

% skip if the 'comb_beat_detectors' option hasn't been specified
if ~sum(strcmp(sig_qual_tools, 'comb_beat_detectors'))
    return
end

%  -- Identify PPG beat detector combinations
beat_detectors = uParams.analysis.beat_detectors;
for detector1_no = 1 : length(beat_detectors)
    for detector2_no = detector1_no : length(beat_detectors)
        % skip if this is the same beat detector
        if detector1_no == detector2_no
            continue
        end
        sig_qual_tools{end+1,1} = ['comb_', beat_detectors{detector1_no}, '_', beat_detectors{detector2_no}];
    end
end

%  -- remove 'comb_beat_detectors'
sig_qual_tools = sig_qual_tools(~strcmp(sig_qual_tools, 'comb_beat_detectors'));

%  -- add in combinations with accel
if sum(strcmp(sig_qual_tools, 'accel'))
    non_accel_tools = sig_qual_tools(~contains(sig_qual_tools, 'accel') & ~strcmp(sig_qual_tools, 'none'));
    for tool_no = 1 : length(non_accel_tools)
        curr_tool = non_accel_tools{tool_no};
        sig_qual_tools{end+1} = ['accel_', curr_tool];
        clear curr_tool
    end
    clear tool_no
end
clear non_accel_tools

end                

function uParams = setup_default_options(uParams, options)
% Setups up the options for the analysis, using the values provided by the 
% user, and default values for options not provided by the user.

%% Specify all the options
option_names = {'beat_detectors'; 'do_downsample'; 'downsample_freq'; 'redo_analysis'};

%% Specify the setting for each option

% setup structure in which to store the settings for the analysis
uParams.analysis = struct;

% cycle through each option in turn
for s = 1 : length(option_names)
    % - check to see whether a setting has been provided
    if sum(strcmp(fieldnames(options),option_names{s}))
        % if so, then use the user-specified option
        eval(['uParams.analysis.' option_names{s} ' = options.' option_names{s} ';']);
    else
        % otherwise, use default option
        switch option_names{s}
            case 'beat_detectors'
                default_setting = {'IMS', 'AMPD', 'MSPTD', 'ABD', 'qppg', 'HeartPy', 'COppg', 'Pulses'};
            case 'do_downsample'
                default_setting = 0;
            case 'downsample_freq'
                default_setting = nan;
            case 'redo_analysis'
                default_setting = 0;
        end
        % store this setting in the uParams.analysis structure:
        eval(['uParams.analysis.' option_names{s} ' = default_setting;']);
        clear default_setting
    end
end
clear s

%% Check that options are compatible
if uParams.analysis.do_downsample && isnan(uParams.analysis.downsample_freq)
    error('\n PPG downsampling has been enabled,\n but no downsampling frequency has been specified.\n Please specify a downsampling frequency\n as explained in the header.');
end

end

function display_startup_message
% Displays a starup message, including details of the licence

licence_details = ['\n\n PPG_BEAT_DETECTOR_ASSESSMENT  Copyright (C) 2022  Peter H. Charlton',...
    '\n This program comes with ABSOLUTELY NO WARRANTY', ... 
    '\n and is available under the GNU public license.', ...
    '\n For details see the accompanying LICENSE.txt file.\n\n'];

fprintf(licence_details)

end

function data = load_data(uParams)
% Load raw dataset from Matlab file for analysis

% Load data from file
data = load_data_from_file(uParams);

% Prepare data for analysis
data = prepare_data_for_analysis(data, uParams);

end

function data = load_data_from_file(uParams)
% Load collated dataset from file

% Check to see whether the specified file exists
if ~exist(uParams.paths.dataset_file, 'file')
    error('\nCouldn''t find data file.\nPlease check that the specified path is correct:\n%s\nIf not, then please either\n (i) provide the filepath as an input in ''options.dataset_file'', or\n (ii) write the correct filepath in the ''setup_up'' function in this script.', uParams.paths.dataset_file);
end

% load data
fprintf('\n - Loading data from: %s', uParams.paths.dataset_file)
load(uParams.paths.dataset_file);

end

function data = prepare_data_for_analysis(data, uParams)
% Prepare the dataset for analysis

% - rename PPG in Vortal
if contains(uParams.dataset, 'vortal')
    for s = 1 : length(data)
        data(s).ppg = data(s).ppgfraw;
    end
end
% - rename ECG in Capnobase and mimic
if contains(uParams.dataset, 'capnobase') || contains(uParams.dataset, 'mimic') || contains(uParams.dataset, 'bidmc')
    for s = 1 : length(data)
        data(s).ecg = data(s).ekg;
        if contains(uParams.dataset, 'capnobase')
            data(s).ecg.rpeaks = data(s).ecg.pk;
        end
    end
end
% - reduce duration of signals if necessary
for s = 1 : length(data)
    ecg_durn = length(data(s).ecg.v)/ data(s).ecg.fs;
    ppg_durn = length(data(s).ppg.v)/ data(s).ppg.fs;
    if ecg_durn > ppg_durn
        final_el = round(ppg_durn*data(s).ecg.fs);
        data(s).ecg.v = data(s).ecg.v(1:final_el);
    end
    if ppg_durn > ecg_durn
        final_el = round(ecg_durn*data(s).ppg.fs);
        data(s).ppg.v = data(s).ppg.v(1:final_el);
    end
end
clear ecg_durn s final_el ppg_durn

% % - remove AF or non-AF data if necessary
% if contains(uParams.dataset, 'mimic')
%     do_af = 0;
%     subjs_to_keep = false(length(data),1);
%     for s = 1 : length(data)
%         if do_af
%             subjs_to_keep(s) = data(s).fix.af_status;
%         else
%             subjs_to_keep(s) = ~data(s).fix.af_status;
%         end
%     end
%     data = data(subjs_to_keep);
% end

% - reduce duration of both ECG and PPG signals to 10 min if necessary
if contains(uParams.dataset, 'mimic')
    durn = 10*60;
    for s = 1 : length(data)
        no_ecg_els = data(s).ecg.fs*durn+1;
        data(s).ecg.v = data(s).ecg.v(1:no_ecg_els);
        no_ppg_els = data(s).ppg.fs*durn+1;
        data(s).ppg.v = data(s).ppg.v(1:no_ppg_els);
    end
    clear s
end

% temporary
%data = data(1:5);

end

function detect_beats_in_ppg_signals(uParams)

fprintf('\n--- Detecting beats in PPG signals')

% cycle through each subject
for subj_no = 1 : uParams.dataset_details.no_subjs
    fprintf('\n - Identifying beats in subj %d: ', subj_no);
    
    %% See if there are any remaining beat detectors to be applied
    % - create filepath
    file_type = 'ppg_beats';
    filepath = create_proc_filepath(uParams, subj_no, file_type);
    % - load file
    if exist(filepath, 'file') && ~uParams.analysis.redo_analysis
        load(filepath);  % loads 'ppg_beats_inds' variable
        beat_detectors_used = fieldnames(ppg_beats_inds);
        beat_detectors_remaining = setxor(beat_detectors_used, uParams.analysis.beat_detectors);
        % - skip this subject if there are no more beat detectors to use
        if isempty(beat_detectors_remaining)
            fprintf('all done');
            continue
        end
    else
        beat_detectors_remaining = uParams.analysis.beat_detectors;
        ppg_beats_inds = struct;
    end
    
    %% Load raw dataset from Matlab file for analysis
    if ~exist('data', 'var')
        data = load_data(uParams);
    end
    
    %% Extract PPG data for this subject
    S = data(subj_no).ppg;
    
    %% Identify periods of no signal
    S = identify_periods_of_no_signal(S, uParams);
    
    %% Pre-process PPG signal
    S = preprocess_ppg_signal(S, uParams);
    
    %% Detect beats in PPG signal
    [ppg_beats_inds, ppg_exc_log] = detect_beats_in_ppg_signal(S, beat_detectors_remaining, ppg_beats_inds, uParams);
    % To check:    plot(S.v), hold on, pks = ppg_beats_inds.IMS; plot(pks, S.v(pks), 'or')
    clear S
    
    %% Save indices of detected beats
    save(filepath, 'ppg_beats_inds', 'ppg_exc_log');
    clear filepath ppg_beats_inds beat_detectors_remaining ppg_exc_log
    
end
clear subj_no

end

function filepath = create_proc_filepath(uParams, subj_no, file_type)

if ~isnan(subj_no)
    eval(sprintf('filepath = ''%s%s_%s.mat'';', uParams.paths.processing_folder, num2str(subj_no,'%04.f'), file_type));
else
    eval(sprintf('filepath = ''%s%s.mat'';', uParams.paths.processing_folder, file_type));
end

end

function S = preprocess_ppg_signal(S, uParams)
% Downsample and filter PPG signal prior to beat detection

% - setup PulseAnalyse options for downsampling
if uParams.analysis.do_downsample
    options.downsample = 1;
    options.downsample_freq = uParams.analysis.downsample_freq;
else
    options.downsample = 0;
end

% Downsample 'no_signal'
if options.downsample
    if rem(S.fs, options.downsample_freq) == 0
        % Downsample
        ds_factor = S.fs/options.downsample_freq;
        S.no_signal = downsample(S.no_signal, ds_factor);
    else
        % Resample
        S.no_signal = resample(double(S.no_signal), options.downsample_freq, S.fs);
        S.no_signal = logical(S.no_signal);
    end
end

% - setup PulseAnalyse options for filtering
options.elim_high_freqs = uParams.analysis.filtering.elim_high_freqs;
options.elim_low_freqs = uParams.analysis.filtering.elim_low_freqs;
options.calc_pw_inds = 0;
options.do_plot = 0;
options.do_beats = 0;
options.do_quality = 0;

% Downsample and filter PPG signal
[~, ~, ~, sigs] = PulseAnalyse(S, options);
S.v = sigs.filt;
S.fs = sigs.fs;
clear options

end

function S = identify_periods_of_no_signal(S, uParams)
% identifies periods when there is no signal (either flat line or nans)

% based on: https://www.mathworks.com/matlabcentral/answers/382011-how-to-count-the-number-of-consecutive-identical-elements-in-both-the-directions-in-a-binary-vecto

% - identify elements which are not nans
not_nan_log = [1;~isnan(S.v(1:end-1));1];
% - identify elements which are not on a flat line
not_flat_line_log = [1;diff(S.v)~=0;1];
% - refine the flat line elements to only include those on a flat line at the max or min envelope of the signal
env_t = 1; % in secs
env_samps = round(env_t*S.fs);
upper_env = movmax(S.v, env_samps);
lower_env = movmin(S.v, env_samps);
on_env = S.v == upper_env | S.v == lower_env;
not_flat_line_log(~not_flat_line_log & [0; ~on_env]) = 1;
% - identify periods of no signal
periods_of_no_signal = diff(find(not_flat_line_log & not_nan_log));
% - only retain periods of no signal which last longer than the threshold duration
durn_of_periods_of_no_signal = repelem(periods_of_no_signal, periods_of_no_signal);
S.no_signal = durn_of_periods_of_no_signal > round(S.fs*uParams.analysis.durn_flat_line);

end

function [ppg_beats_inds, ppg_exc_log] = detect_beats_in_ppg_signal(S, beat_detectors_remaining, ppg_beats_inds, uParams)

% set up variable to note any periods of PPG to be excluded from the analysis
ppg_exc_log = false(size(S.v));

% identify beats in PPG
for beat_detector_no = 1 : length(beat_detectors_remaining)
    
    curr_beat_detector = beat_detectors_remaining{beat_detector_no};    
    fprintf('%s, ', curr_beat_detector);
    
    % - setup PulseAnalyse options for beat detection
    options.calc_pw_inds = 0;
    options.do_plot = 0;
    options.do_beat_filter = 0;
    options.do_filter = 0;
    options.do_quality = 0;
    options.beat_detector = curr_beat_detector;
    
    % - detect beats in this subject's PPG signals using this peak detector
    if ~uParams.analysis.do_wins
        % Without windowing
        [~, ~, pulses, ~] = PulseAnalyse(S, options);
        if strcmp(uParams.analysis.ppg_fid_pt, 'pk')
            eval(['ppg_beats_inds.' curr_beat_detector ' = pulses.peaks;']);
        elseif strcmp(uParams.analysis.ppg_fid_pt, 'onset')
            eval(['ppg_beats_inds.' curr_beat_detector ' = pulses.onsets;']);
        elseif strcmp(uParams.analysis.ppg_fid_pt, 'mid_pt')
            eval(['ppg_beats_inds.' curr_beat_detector ' = pulses.mid_amps;']);
        end
    else
        % With windowing
        win_starts = 0:(uParams.analysis.win_durn-uParams.analysis.win_overlap):(((length(S.v)-1)/S.fs)-uParams.analysis.win_durn);
        t = 0:(1/S.fs):((length(S.v)-1)/S.fs);
        temp_beats = [];
        for win_no = 1 : length(win_starts)
            
            % identify elements corresponding to this window
            curr_start = win_starts(win_no);
            curr_end = curr_start+uParams.analysis.win_durn;
            rel_els = t>= curr_start & t<= curr_end;
            
            % skip this window if any of it didn't contain a signal
            if sum(S.no_signal(rel_els))>0
                ppg_exc_log(rel_els) = true;
                continue
            end
            
            % If this is the first window then exclude the first part from the analysis
            if win_no == 1
                els_to_exc = t>= curr_start & t< (curr_start+uParams.analysis.win_overlap);
                ppg_exc_log(els_to_exc) = true;
            end
            
            % identify beats in this window
            rel_S.v = S.v(rel_els);
            rel_S.fs = S.fs;
            [~, ~, pulses, ~] = PulseAnalyse(rel_S, options);
            
            % store beats
            if strcmp(uParams.analysis.ppg_fid_pt, 'pk')
                temp_beats = [temp_beats; pulses.peaks+find(rel_els,1)-1];
            elseif strcmp(uParams.analysis.ppg_fid_pt, 'onset')
                temp_beats = [temp_beats; pulses.onsets+find(rel_els,1)-1];
            elseif strcmp(uParams.analysis.ppg_fid_pt, 'mid_pt')
                temp_beats = [temp_beats; pulses.mid_amps+find(rel_els,1)-1];
            end
            clear pulses curr_start curr_end rel_els rel_S
        end
        clear win_no curr_start curr_end rel_els
        
        % remove repeated beat detections
        temp_beats = sort(temp_beats);
        if ~isempty(temp_beats)
            repeated_beats = [0; diff(temp_beats)< round(uParams.analysis.tol_window*S.fs)];
            temp_beats = temp_beats(~repeated_beats);
            clear repeated_beats
        end
        
        % add in last bit which wasn't captured in the windows to the exclusion logical:
        rel_els = t>= (win_starts(end) + uParams.analysis.win_durn);
        ppg_exc_log(rel_els) = true;  
        clear win_starts t rel_els
        
        % store beats
        eval(['ppg_beats_inds.' curr_beat_detector ' = temp_beats;']);
        clear temp_beats win_no t win_starts
    end
    clear pulses do_wins options curr_beat_detector
    
end
clear beat_detector_no

end

function assess_quality_of_ppg_signals(uParams)

fprintf('\n--- Assessing quality of PPG signals')

% cycle through each subject
for subj_no = 1 : uParams.dataset_details.no_subjs
    fprintf('\n - Assessing quality for subj %d: ', subj_no);
    
    %% See if there are any remaining quality tools to be applied
    % - create filepath
    file_type = 'ppg_qual';
    filepath = create_proc_filepath(uParams, subj_no, file_type);
    
    if exist(filepath, 'file') && ~uParams.analysis.redo_analysis
        % - if processing has been started before, then find out which tools haven't yet been applied to this subject's data
        load(filepath);  % loads 'ppg_qual' variable
        qual_tools_used = fieldnames(ppg_qual);
        qual_tools_remaining = setxor(qual_tools_used, uParams.analysis.sig_qual_tools);
        % - skip this subject if there are no more quality assessment tools to use
        if isempty(qual_tools_remaining)
            fprintf('all done');
            continue
        end
    else
        % - otherwise, use all the tools
        qual_tools_remaining = uParams.analysis.sig_qual_tools;
        ppg_qual = struct;
    end
    
    % Move accel-based tools to the end (because some of these rely on non-accel based tools)
    accel_based_tools = contains(qual_tools_remaining, 'accel');
    qual_tools_remaining = qual_tools_remaining([find(~accel_based_tools); find(accel_based_tools)]);
    clear accel_based_tools
    
    % cycle through each tool
    for tool_no = 1 : length(qual_tools_remaining)
        curr_tool = qual_tools_remaining{tool_no};
    
        %% Load raw dataset from Matlab file for analysis
        if ~exist('data', 'var')
            data = load_data(uParams);
        end
        
        %% use this tool to assess quality
        
        if strcmp(curr_tool, 'none')
            % - if no tool, then just label all samples as of high quality.
            curr_qual.v = true(length(data(subj_no).ppg.v),1);
            curr_qual.fs = data(subj_no).ppg.fs;
            eval(['ppg_qual.' curr_tool ' = curr_qual;']);
            
        elseif strcmp(curr_tool, 'accel')
            % - use accelerometry signal at same site as PPG measurement
            fprintf('accel, ');
            curr_qual = assess_qual_using_accel(data(subj_no).acc_ppg_site, data(subj_no).ppg, subj_no, curr_tool, uParams);
            
            % store results
            eval(['ppg_qual.' curr_tool ' = curr_qual;']);
            clear curr_qual curr_tool
            
        elseif strcmp(curr_tool(1:5), 'comb_')
            % - use combinations of PPG beat detectors
            
            %  -- Load PPG beat detections
            file_type = 'ppg_beats';
            loadpath = create_proc_filepath(uParams, subj_no, file_type);
            load(loadpath); % load 'ppg_beat_inds' variable
            
            %  -- Assess qual using each combination of PPG beat detectors
            temp = strfind(curr_tool, '_');
            detector1 = curr_tool(temp(1)+1:temp(2)-1);
            detector2 = curr_tool(temp(2)+1:end);
            for detector_no = 1 : 2
                eval(['detector' num2str(detector_no) '_ppg_beats_inds = ppg_beats_inds.' eval(['detector', num2str(detector_no)]), ';']);
            end
            fprintf([curr_tool ', ']);
            curr_qual = assess_qual_using_comb_beat_detectors(data(subj_no).ppg, detector1_ppg_beats_inds, detector2_ppg_beats_inds, uParams);
            clear detector1 detector2 detector1_ppg_beats_inds detector2_ppg_beats_inds
            
            % store results
            eval(['ppg_qual.' curr_tool ' = curr_qual;']);
            clear curr_qual curr_tool
            
        elseif strcmp(curr_tool(1:6), 'accel_')
            % use accel with another tool
            % - identify the other tool
            temp = strfind(curr_tool, '_');
            other_tool = curr_tool(temp(1)+1:end);
            accel_qual = ppg_qual.accel;
            eval(['other_qual = ppg_qual.' other_tool ';']);
            if ~isequal(other_qual.fs, accel_qual.fs)
                error('\n These sampling freqs should be the same')
            end
            curr_qual.fs = other_qual.fs;
            curr_qual.v = other_qual.v & accel_qual.v;
            
            % store results
            eval(['ppg_qual.' curr_tool ' = curr_qual;']);
            clear curr_qual curr_tool accel_qual other_qual temp
        end
        
    end
    clear tool_no
    
    %% Save quality results for this subject
    save(filepath, 'ppg_qual');
    clear filepath ppg_qual qual_tools_remaining
    
end
clear subj_no

end

function curr_qual = assess_qual_using_accel(acc_ppg_site, ppg, subj_no, curr_tool, uParams)

%% Obtain mean absolute deviation (MAD) signal
%  - create time vector
acc_ppg_site.t = [0:length(acc_ppg_site.v)-1]./acc_ppg_site.fs; acc_ppg_site.t = acc_ppg_site.t(:);
%  - segment into 5 sec windows
win_durn = 5; % 5 secs as specified in ﻿https://doi.org/10.1371/journal.pone.0164045 and ﻿https://doi.org/10.1111/cpf.12127
win_starts = 0:win_durn:acc_ppg_site.t(end)-win_durn;
mad_sig.v = nan(length(win_starts),1);
mad_sig.fs = 1/5; % one sample every 5 secs
mad_sig.t = [0:length(mad_sig.v)-1]./mad_sig.fs; mad_sig.t = mad_sig.t(:);
%  - calculate MAD value for each 5 sec window
for win_no = 1 : length(win_starts)
    rel_els = acc_ppg_site.t>=win_starts(win_no) & acc_ppg_site.t<(win_starts(win_no)+win_durn);
    deviations = acc_ppg_site.v(rel_els)-mean(acc_ppg_site.v(rel_els));
    mad_sig.v(win_no) = mean(abs(deviations));
end
clear win_no deviations rel_els win_starts win_durn

%% Assess quality based on MAD signal
mad_sig.hq = mad_sig.v< 16.7; % threshold from p.68 of ﻿https://doi.org/10.1111/cpf.12127

%% Create vector of qualities at PPG's original sampling freq
curr_qual = create_vector_of_quals_at_ppg_freq(ppg, mad_sig);

end

function curr_qual = assess_qual_using_comb_beat_detectors(ppg, detector1_ppg_beats_inds, detector2_ppg_beats_inds, uParams)

%% Convert beat indices to timings
for detector_no = 1 :2
    eval(['detector' num2str(detector_no) '_t = [detector' num2str(detector_no) '_ppg_beats_inds-1]./ppg.fs;']);
end
clear detector1_ppg_beats_inds detector2_ppg_beats_inds detector_no

%% Select order of detectors (according to assumed sensitivity and PPV)
if length(detector1_t) > length(detector2_t)
    sens_detector_t = detector1_t;
    high_ppv_detector_t = detector2_t;
else
    sens_detector_t = detector2_t;
    high_ppv_detector_t = detector1_t;
end
clear detector1_t detector2_t

%% find differences in timings
diff_matrix = repmat(sens_detector_t, [1, length(high_ppv_detector_t)]) - high_ppv_detector_t';
% - Find minimum differences
[min_abs_diff, used_detector1_els] = min(abs(diff_matrix));
min_abs_diff = min_abs_diff(:);
% - Identify correctly identified beats
beat_correct_log = min_abs_diff<uParams.analysis.tol_window;
clear min_abs_diff diff_matrix
% - identify times of incorrect detections
incorrect_high_ppv_detector_t = high_ppv_detector_t(~beat_correct_log);
used_detector1_log = false(length(sens_detector_t),1);
used_detector1_log(used_detector1_els) = true;
incorrect_sens_detector_t = sens_detector_t(~used_detector1_log);
incorrect_t = unique([incorrect_sens_detector_t; incorrect_high_ppv_detector_t]);
clear incorrect_sens_detector_t incorrect_high_ppv_detector_t used_detector1_log used_detector1_els
% - List agreed beat detections
agreed_beats.t = high_ppv_detector_t(beat_correct_log);

%% assess quality of each beat detection

% - High qual if: (i) previous IBI < 2 s; (ii) next IBI < 2 s; 
curr_IBI = [inf; diff(agreed_beats.t); inf];
high_qual_IBI_log = curr_IBI(1:end-1) < 2 & curr_IBI(2:end) < 2;
clear curr_IBI

% - High qual if: (iii) no incorrect_t in previous 5 s
if ~isempty(incorrect_t)
    diff_matrix = agreed_beats.t' - repmat(incorrect_t, [1, length(agreed_beats.t)]);
    diff_matrix(diff_matrix<0) = inf;
    min_diff = min(diff_matrix);
    min_diff = min_diff(:);
    clear diff_matrix
else
    min_diff = zeros(length(beat_correct_log),1);
end
clear incorrect_t beat_correct_log

% - Identify high quality beats (based on IBIs and whether or not there was an incorrect detection in the previous 5s)
agreed_beats.hq = high_qual_IBI_log & min_diff>uParams.analysis.window_durn_for_incorrect_beats;
clear min_diff high_qual_IBI_log

%% Create vector of qualities at PPG's original sampling freq
curr_qual = create_vector_of_quals_at_ppg_freq(ppg, agreed_beats);

end

function curr_qual = create_vector_of_quals_at_ppg_freq(ppg, qual_vector)

% - create time vector for ppg

% - say all low quality if there are no entries in the qual vector
if isempty(qual_vector.t)
    curr_qual.fs = ppg.fs;
    curr_qual.v = false(length(ppg.v),1);
    return
end

% - resample at original sampling freq
ppg.t = [0:length(ppg.v)-1]./ppg.fs; ppg.t = ppg.t(:);
curr_qual.v = interp1(qual_vector.t, double(qual_vector.hq), ppg.t, 'previous');
curr_qual.fs = ppg.fs;
% - convert to logical
curr_qual.v(isnan(curr_qual.v)) = 0;
curr_qual.v = logical(curr_qual.v);

end

function detect_beats_in_ecg_signals(uParams)

fprintf('\n--- Detecting beats in ECG signals')

% cycle through each subject
for subj_no = 1 : uParams.dataset_details.no_subjs
    fprintf('\n - Identifying ECG beats in subj %d: ', subj_no);
    
    %% See if there are any remaining beat detectors to be applied
    % - create filepath
    file_type = 'ecg_beats';
    filepath = create_proc_filepath(uParams, subj_no, file_type);
    % - skip if this file exists (and therefore the processing has been done)
    if exist(filepath, 'file') && ~uParams.analysis.redo_analysis
        fprintf('all done');
        continue
    end
    
    %% Load raw dataset from Matlab file for analysis
    if ~exist('data', 'var')
        data = load_data(uParams);
    end
    
    %% Extract ECG data for this subject
    S = data(subj_no).ecg;
    
    %% Identify periods of no signal
    S = identify_periods_of_no_signal(S, uParams);
    
    %% Detect beats in ECG signal
    [ecg_beats_inds, ecg_exc_log] = detect_beats_in_ecg_signal(S, uParams);
    % To check:         plot(S.v), hold on, pks = ecg_beats_inds; plot(pks, S.v(pks), 'or')
    clear S
    
    %% Save indices of detected beats
    save(filepath, 'ecg_beats_inds', 'ecg_exc_log');
    clear filepath ecg_beats_inds ecg_exc_log
    
end
clear subj_no

end

function [ecg_beats_inds, ecg_exc_log] = detect_beats_in_ecg_signal(S, uParams)
% detects beats in an ECG signal

ecg_exc_log = false(size(S.v));

%% check to see whether there are already annotations of R-peaks in the dataset
if sum(strcmp(fieldnames(S), 'rpeaks'))
    ecg_beats_inds = S.rpeaks(:);
    ecg_beats_inds = unique(ecg_beats_inds); % inserted to eliminate any repeated annotations (because I found one once)
    return
end

%% if not, then identify R-peaks

do_orig = false;

if do_orig
    if ~uParams.analysis.do_wins
        % Without windowing
        [~, ecg_beats_inds, ~]  = rpeakdetect_pc(S.v,S.fs);
    else
        % With windowing
        win_starts = 0:(uParams.analysis.win_durn-uParams.analysis.win_overlap):(((length(S.v)-1)/S.fs)-uParams.analysis.win_durn);
        t = 0:(1/S.fs):((length(S.v)-1)/S.fs);
        ecg_beats_inds = [];
        for win_no = 1 : length(win_starts)
            
            % identify elements corresponding to this window
            curr_start = win_starts(win_no);
            curr_end = curr_start+uParams.analysis.win_durn;
            rel_els = t>= curr_start & t<= curr_end;
            
            % skip this window if any of it didn't contain a signal
            if sum(S.no_signal(rel_els))
                ecg_exc_log(rel_els) = true;
                continue
            end
            
            % identify beats in this window
            rel_S.v = S.v(rel_els);
            rel_S.fs = S.fs;
            [~, curr_ecg_beat_inds, ~]  = rpeakdetect_pc(rel_S.v,rel_S.fs);
            
            % store beats
            ecg_beats_inds = [ecg_beats_inds; curr_ecg_beat_inds+find(rel_els,1)-1];
            clear curr_ecg_beat_inds curr_start curr_end rel_els rel_S
        end
        clear win_no t win_starts
        
        % remove repeated beat detections
        ecg_beats_inds = sort(ecg_beats_inds);
        repeated_beats = [0; diff(ecg_beats_inds)< round(uParams.analysis.tol_window*S.fs)];
        ecg_beats_inds = ecg_beats_inds(~repeated_beats);
        clear repeated_beats win_no t win_starts
    end
    clear do_wins
else
    
    % Detect heartbeats in an ECG signal, and assess the quality of the beat detections.
    options.verbose = false;
    options.win_durn = uParams.analysis.win_durn;
    options.win_overlap = uParams.analysis.win_overlap;
    options.qrs_tol_window = uParams.analysis.qrs_tol_window;
    [ecg_beats_inds, qual] = ecg_beat_detector(S.v, S.fs, options, S.no_signal);
    ecg_exc_log = 1-qual;
    
%     %% GQRS detector
%     % - write this signal to a WFDB file
%     temp_diffs = diff(sort(S.v));
%     temp_sig = round(S.v/(min(temp_diffs(temp_diffs>0))));
%     tm = (0:(length(S.v)-1))/S.fs;
%     wrsamp(tm(:),temp_sig,'temp',S.fs);
%     
%     % - run gqrs to identify beats, and save to a qrs annotation file
%     gqrs('temp');
%     
%     % - load qrs annotations
%     [ann]=rdann('temp', 'qrs');
%     
%     % - display message about temporary files
%     fprintf('\n temporary files created in current directory')
%     
%     %% JQRS detector
%     HRVparams = InitializeHRVparams(''); HRVparams.Fs = S.fs;
%     QRS = run_qrsdet_by_seg(S.v,HRVparams);
%     QRS = QRS(:);
%     
%     %% Assess quality
%     % - here ann is the reference
%     diff_matrix = repmat(QRS, [1, length(ann)]) - ann';
%     % - Find minimum differences
%     min_abs_diff = min(abs(diff_matrix), [], 2);
%     % - Identify correctly identified beats
%     beat_correct_log = min_abs_diff<(uParams.analysis.qrs_tol_window*S.fs);
%     ecg_beats_inds = QRS(beat_correct_log);
%     % remove any repeated beat detections
%     ecg_beats_inds = sort(ecg_beats_inds);
%     repeated_beats = [0; diff(ecg_beats_inds)< round(uParams.analysis.tol_window*S.fs)];
%     ecg_beats_inds = ecg_beats_inds(~repeated_beats);
%     clear repeated_beats
%         
%     if uParams.analysis.do_wins
%         
%         % Perform quality assessment with windowing
%         win_starts = 0:(uParams.analysis.win_durn-uParams.analysis.win_overlap):(((length(S.v)-1)/S.fs)-uParams.analysis.win_durn);
%         t = 0:(1/S.fs):((length(S.v)-1)/S.fs);
%         t_disagree_beats = t(QRS(~beat_correct_log));
%         for win_no = 1 : length(win_starts)
%             
%             % identify elements corresponding to this window
%             curr_start = win_starts(win_no);
%             curr_end = curr_start+uParams.analysis.win_durn;
%             
%             % skip this window if any of it didn't contain a signal
%             rel_els = t>= curr_start & t<= curr_end;
%             if sum(S.no_signal(rel_els))
%                 ecg_exc_log(rel_els) = true;
%                 continue
%             end
%             
%             % skip this window if it contained any disagreements
%             curr_win_disagree_inds = t_disagree_beats >= curr_start & t_disagree_beats <= curr_end;
%             if sum(curr_win_disagree_inds)
%                 ecg_exc_log(rel_els) = true;
%                 continue
%             end
%             
%             clear curr_start curr_end rel_els curr_win_disagree_inds
%         end
%         clear win_no t win_starts t_disagree_beats
%         
%     end
end

end

function time_align_ppg_beats(uParams)

fprintf('\n--- Aligning times of PPG and ECG beats')

% cycle through each subject
for subj_no = 1 : uParams.dataset_details.no_subjs
    fprintf('\n - Aligning beats for subj %d: ', subj_no);
    
    %% See if this has been done already
    % - create filepath
    file_type = 'ecg_beats_aligned';
    filepath = create_proc_filepath(uParams, subj_no, file_type);
    % - load file
    if exist(filepath, 'file') && ~uParams.analysis.redo_analysis
        load(filepath);  % loads 'ecg_beats_a_inds' variable
        beat_detectors_used = fieldnames(ecg_beats_a_inds);
        beat_detectors_remaining = setxor(beat_detectors_used, uParams.analysis.beat_detectors);
        % - skip this subject if there are no more beat detectors to use
        if isempty(beat_detectors_remaining)
            fprintf('all done');
            continue
        end
    else
        beat_detectors_remaining = uParams.analysis.beat_detectors;
        ecg_beats_a_inds = struct;
    end
        
    %% Load PPG Beats and ECG Beats files for analysis
    file_type = 'ppg_beats';
    loadpath = create_proc_filepath(uParams, subj_no, file_type);
    load(loadpath)
    file_type = 'ecg_beats';
    loadpath = create_proc_filepath(uParams, subj_no, file_type);
    load(loadpath)
    clear loadpath file_type
    
    %% Find time-aligned PPG beat timings for each PPG beat detector
    for detector_no = 1 : length(beat_detectors_remaining)
        curr_detector = beat_detectors_remaining{detector_no};
        fprintf([curr_detector, ', ']); 
    
        %% Time-align ECG with PPG beats obtained using this beat detector
        eval(['curr_ppg_beats_inds = ppg_beats_inds.' curr_detector ';']);
        [curr_ecg_beats_a_inds, ecg_exc_a_log] = time_align_ppg_beats_with_ecg_beats(curr_ppg_beats_inds, ecg_beats_inds, ecg_exc_log, uParams);
        
        % store results
        eval(['ecg_beats_a_inds.' curr_detector ' = curr_ecg_beats_a_inds;']);
        clear curr_ecg_beats_a_inds curr_detector curr_ecg_beats_inds
        
    end
    clear detector_no ppg_beats_inds ecg_beats_inds
    
    %% Save quality results for this subject
    save(filepath, 'ecg_beats_a_inds', 'ecg_exc_a_log');
    clear filepath ecg_beats_a_inds beat_detectors_remaining ecg_exc_a_log
    
end
clear subj_no

end

function [ecg_beats_a_inds, ecg_exc_a_log] = time_align_ppg_beats_with_ecg_beats(ppg_beats_inds, ecg_beats_inds, ecg_exc_log, uParams)

% in the future this could be adjusted to include exclusion of low quality ECG beats with ecg_exc_log.

% - skip if there are no detected PPG beats
if isempty(ppg_beats_inds)
    ecg_beats_a_inds = ecg_beats_inds;
    ecg_exc_a_log = ecg_exc_log;
    return
end

% Transform beat indices to timings
ecg_beats_t = [ecg_beats_inds-1]./uParams.dataset_details.ecg_fs;
ppg_beats_t = [ppg_beats_inds-1]./uParams.dataset_details.ppg_fs_ds;

% - Calculate time differences between each pair of PPG and ECG beats
diff_matrix = repmat(ppg_beats_t, [1, length(ecg_beats_t)]) - ecg_beats_t';

% - calculate performance for a range of lags
lag_ts = (-1*uParams.analysis.max_lag):uParams.analysis.lag_int:uParams.analysis.max_lag;
curr_perf = nan(size(lag_ts));
for lag_t_el = 1:length(lag_ts)
    lag_t = lag_ts(lag_t_el);
    % - Find minimum differences for this lag
    min_abs_diff = min(abs(diff_matrix-lag_t));
    min_abs_diff = min_abs_diff(:);
    % - Identify correctly identified beats when using this lag
    beat_correct_log = min_abs_diff<uParams.analysis.tol_window;
    curr_perf(lag_t_el) = sum(beat_correct_log);
end
clear min_abs_diff beat_correct_log lag_t lag_t_el

% - identify the lag which gives the best performance
temp = max(curr_perf);
rel_lag_t_els = find(curr_perf == temp);  % because more than one lag_t can result in this f1
[~, temp2] = min(abs(lag_ts(rel_lag_t_els)));
rel_lag_t = lag_ts(rel_lag_t_els(temp2));
rel_lag_ecg_samps = round(rel_lag_t*uParams.dataset_details.ecg_fs);
clear temp2 temp curr_perf rel_lag_t_els lag_ts

% - calculate the aligned indices of PPG beats
ecg_beats_a_inds = ecg_beats_inds+rel_lag_ecg_samps;

% - calculate time-shifted ecg exclusion log
if rel_lag_ecg_samps > 0
    ecg_exc_a_log = [true(rel_lag_ecg_samps,1); ecg_exc_log(1:end-rel_lag_ecg_samps)];
elseif rel_lag_ecg_samps < 0
    ecg_exc_a_log = [ecg_exc_log(abs(rel_lag_ecg_samps)+1:end); true(abs(rel_lag_ecg_samps),1)];
else
    ecg_exc_a_log = ecg_exc_log;
end

% To check:
% plot(ppg_beats_t, zeros(size(ppg_beats_t)), '*k'), hold on, plot(ecg_beats_t, zeros(size(ecg_beats_t)), 'or'), plot(ppg_beats_t, ones(size(ppg_beats_t)), '*k'), hold on, plot(ecg_beats_t+rel_lag_t, ones(size(ecg_beats_t)), 'or'), ylim([-1, 2])

end

function ints_ts = generate_beat_signal(beat_timings, up)

ints.t = beat_timings(1:end-1);
ints.v = diff(beat_timings);
els_to_exclude = isnan(ints.t) | isnan(ints.v);
ints.t = ints.t(~els_to_exclude);
ints.v = ints.v(~els_to_exclude);

ints_ts.t = 0:(1/up.analysis.interpolation_fs):max(beat_timings);
ints_ts.v = interp1(ints.t, ints.v, ints_ts.t, 'previous');

ints_ts.t = round(ints_ts.t*100)/100; % assumes a max sampling freq of 100 Hz

ints_ts.t = ints_ts.t(:);
ints_ts.v = ints_ts.v(:);
end

function hr_ts = generate_hr_signal(beat_timings, beat_hrs, signal_start_t, signal_end_t, up)

hrs.t = beat_timings;
hrs.v = beat_hrs;
els_to_exclude = isnan(hrs.t) | isnan(hrs.v);
hrs.t = hrs.t(~els_to_exclude);
hrs.v = hrs.v(~els_to_exclude);

hr_ts.t = signal_start_t:(1/up.analysis.interpolation_fs):signal_end_t;

% skip if there are less than two HRs to interpolate (the minimum number required by this function)
if length(hrs.t)< 2
    hr_ts.v = nan(size(hr_ts.t));
else
    hr_ts.v = interp1(hrs.t, hrs.v, hr_ts.t, 'previous');
end

hr_ts.t = round(hr_ts.t*100)/100; % assumes a max sampling freq of 100 Hz

% make into column vectors
hr_ts.t = hr_ts.t(:);
hr_ts.v = hr_ts.v(:);

end

function qual_ts = generate_qual_signal(beat_timings, beat_qual, up)

qual.t = beat_timings;
qual.v = double(beat_qual);
%els_to_exclude = isnan(ints.t) | isnan(ints.v);
%ints.t = ints.t(~els_to_exclude);
%ints.v = ints.v(~els_to_exclude);

qual_ts.t = 0:(1/up.analysis.interpolation_fs):max(beat_timings);
qual_ts.v = interp1(qual.t, qual.v, qual_ts.t, 'previous');

qual_ts.t = round(qual_ts.t*100)/100; % assumes a max sampling freq of 100 Hz

% make into column vectors
qual_ts.t = qual_ts.t(:);
qual_ts.v = qual_ts.v(:);
end

function [R_amp, R_index, S_amp]  = rpeakdetect_pc(data,samp_freq,thresh)

% [R_amp, R_index, S_amp]  = rpeakdetect(data, samp_freq, thresh); 
% R_amp == amplitude of R peak in bpf data & S_amp == amplitude of following minmum.
% samp_freq == sampling frequency in Hz 
% The 'triggering' threshold 'thresh' for the peaks in the 'integrated'  
% waveform is 0.2 by default.
%
% A batch QRS detector based upon that of Pan, Hamilton and Tompkins:
% J. Pan \& W. Tompkins - A real-time QRS detection algorithm 
% IEEE Transactions on Biomedical Engineering, vol. BME-32 NO. 3. 1985.
% P. Hamilton \& W. Tompkins. Quantitative Investigation of QRS 
% Detection  Rules Using the MIT/BIH Arrythmia Database. 
% IEEE Transactions on Biomedical Engineering, vol. BME-33, NO. 12.1986.
% 
% Similar results reported by the authors above were achieved, without
% having to tune the thresholds on the MIT DB. An online version in C
% has also been written.
%
% Written by G. Clifford gari@ieee.org and made available under the 
% GNU general public license. If you have not received a copy of this 
% license, please download a copy from http://www.gnu.org/
%
% Please distribute (and modify) freely, commenting
% where you have added modifications. 
% The author would appreciate correspondence regarding
% corrections, modifications, improvements etc.
%
% gari@ieee.org
%
% --- 28th Feb 2019, Peter H Charlton ---
% Made the following changes:
%  - Added further comments (an adjusted original ones) to help me understand the code
%  - Split the integration step into two parts: integrating, and taking the
%  median filter of the integrated signal.
%  - Removed time axis and plotting function (including testmode and time outputs)
%  - Removed default sampling frequency
%  - Removed HRV and resp variables
% --- 13th July 2021, Peter H Charlton ---
%  - Made maxval, minval, maxloc, minloc into column vectors
% ---------------------------------------

%% Setup

% --- Threshold ---

%- make threshold default 0.2 -> this was 0.15 on MIT data 
if nargin < 3
   thresh = 0.2;
end

% --- Data Format ---

% - Make into column vector
x = data(:);

% eliminate nans
t = 0:length(x);
x(isnan(x)) = interp1(t(~isnan(x)),x(~isnan(x)),t(isnan(x)));

%% Filter signal

% --- Remove mean ---

x = x-nanmean(x);  % there may still be some nans at the start or end
 
% --- FIR Bandpass filter (if filter exists) ---

% - no filter
bpf=x; % use if there is no filter bank available
% - FS = 128 Hz
if( (samp_freq==128) & (exist('filterECG128Hz')~=0) )
        bpf = filterECG128Hz(x); 
end
% - FS = 256 Hz
if( (samp_freq==256) & (exist('filterECG256Hz')~=0) )
        bpf = filterECG256Hz(x); 
end

%% Process signal

% --- Differentiate ---
dff = diff(bpf);  % now it's one datum shorter than before

% --- Square ---
sqr = dff.*dff;   %
len = length(x)-1; % how long is the new vector now?

% --- Integrate data over window 'd' ---
% - setup window size
d=[1 1 1 1 1 1 1]; % intialise
if (samp_freq>=256) % adapt for higher sampling rates
    d = [ones(1,round(7*samp_freq/256))];
end
% - integrate
integ = filter(d,1,sqr);
% - take median filter over window of 10 samples
mdfint = medfilt1(integ,10);
% - remove filter delay for scanning back through ECG
delay = ceil(length(d)/2);
mdfint = mdfint(delay:length(mdfint));

%% Search

% - Maximum value in middle half of processed signal
max_h = max (mdfint(round(len/4):round(3*len/4)));

% - Identify indices of processed signal which are greater than the threshold proportion of this maximum value
poss_reg = mdfint>(thresh*max_h);

% - Find regions of processed signal which are above this threshold
left  = find(diff([0 poss_reg'])==1); % left boundary, remembering to zero pad at start
right = find(diff([poss_reg' 0])==-1); % right boundary, remembering to zero pad at end

% - identify maximum and minimum values in each region in BPF signal
for(i=1:length(left))
    [maxval(i,1), maxloc(i,1)] = max( bpf(left(i):right(i)) );
    [minval(i,1), minloc(i,1)] = min( bpf(left(i):right(i)) );
    maxloc(i,1) = maxloc(i)-1+left(i); % add offset of present location
    minloc(i,1) = minloc(i)-1+left(i); % add offset of present location
end

% - Assume R-peaks are at maximum values of BPF signal
R_index = maxloc;
R_amp = maxval;

% - Assume S-troughs are at minimum values of BPF signal
S_amp = minval;   %%%% Assuming the S-wave is the lowest amp in the given window

% - check for lead inversion (i.e. do minima precede maxima?)
if (minloc(length(minloc))<maxloc(length(minloc)))      % considering the final beat only
    R_amp = minval;
    S_amp = maxval;
end

end
