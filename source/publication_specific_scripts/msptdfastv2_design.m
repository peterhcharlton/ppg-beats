function msptdfastv2_design
% MSPTDFASTV2_DESIGN  Used for MSPTDfast v.2 design.
%   MSPTDFASTV2_DESIGN analyses the performance
%   (execution time and F1-score) of different configurations of the MSPTD
%   PPG beat detector.
%   
%   # Inputs
%   * ppg_detect_stats.mat - results files from analysis of PPG-DaLiA
%   (all_activities), 
%   lunchbreak data.
%   
%   # Outputs
%   * Plots illustrating the performance of PPG beat detectors.
%   
%   # Exemplary usage
%   
%       msptdfast_cinc_analysis_post_20240701      % runs the analysis
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Author
%   Peter H. Charlton, University of Cambridge, July 2024.
%   
%   # License - MIT
%      Copyright (c) 2024 Peter H. Charlton
%      Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
%      The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
%      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% Setup universal parameters
up = setup_up;

% cycle through each dataset
for dataset_no = 1 : length(up.datasets)

    curr_dataset = up.datasets{dataset_no};

    fprintf('\n ~~~ Generating results plots for %s ~~~', curr_dataset)

    % load data
    eval(['up.paths.stats_file = up.paths.stats_file_' curr_dataset ';']);
    load(up.paths.stats_file)

    % generate plot for each option
    generate_results_plots(res, curr_dataset, up);

end

end

function up = setup_up

up.paths.stats_file_ppg_dalia_all_activities = '/Users/petercharlton/Documents/Data/ppg_dalia/conv_data/proc_data_ppg_dalia_all_activities/ppg_detect_stats.mat';
up.paths.stats_file_mimic_train_all = '/Users/petercharlton/Documents/Data/mimic_perform_train_test_datasets/proc_data_mimic_train_all/ppg_detect_stats.mat';
up.paths.stats_file_capnobase = '/Users/petercharlton/Documents/Data/capnobase_2021/conv_data/proc_data_capnobase/ppg_detect_stats.mat';
up.paths.stats_file_bidmc = '/Users/petercharlton/Documents/Data/BIDMC_dataset/proc_data_bidmc/ppg_detect_stats.mat';
up.paths.plot_folder = '/Users/petercharlton/Library/CloudStorage/GoogleDrive-peterhcharlton@gmail.com/My Drive/Work/Publications/In Preparation/2024 MSPTDfastv2/figures/';

up.datasets = {'ppg_dalia_all_activities', 'mimic_train_all', 'capnobase', 'bidmc'};
up.dataset_names = {'PPG-DaLiA (all activities)', 'MIMIC PERform (training)', 'CapnoBase', 'BIDMC'};

up.options.name = {'lmss'; ...
    'lms_method'; ...
    'lms_scales'; ...
    'fs'; ...
    'win_durn'; ...
    %'algs'; ...
    %'red_cols'; ...
    };
up.options.label = {'LMSs to calculate'; 
    'LMS calculation method'; ...
    'LMS scales used'; ...
    'Sampling Frequency (Hz)'; ...
    'Window Duration (s)'; ...
    %'Algorithm'; ...
    %'LMS columns used'; ...
    };
up.options.values = {{'peaks', 'onsets', 'peaks and onsets'}; ...
    {'vectorised', 'nested loops'}; ...
    {'HRmin=40', 'HRmin=30', 'N/2'}; ...
    {'10', '20', '30', '40', '50', 'orig'}; ...
    {'4', '6', '8', '10', '12'}; ...
    %{'MSPTD', 'MSPTDfastv1.1', 'qppgfast'}; ...
    %{'Reduced set', 'All'}; ...
    };
up.options.type = {'class'; 
    'class'; ...
    'continuous'; ...
    'continuous'; ...
    'continuous'; ...
    %'class'; ...
    %'class'; ...
    };
up.options.alg_abbr = {{'MSPTDPC1', 'MSPTDPC2', 'MSPTDPC5'}; ...
    {'MSPTDPC3', 'MSPTDPC5'}; ...
    {'MSPTDPC12', 'MSPTDPC4', 'MSPTDPC5'}; ...
    {'MSPTDPC8', 'MSPTDPC7', 'MSPTDPC6', 'MSPTDPC15', 'MSPTDPC16', 'MSPTDPC5'}; ...
    {'MSPTDPC9', 'MSPTDPC10', 'MSPTDPC5', 'MSPTDPC11', 'MSPTDPC14'}; ...
    %{'MSPTD', 'MSPTDfastv1', 'qppgfast'}; ...
    %{'MSPTDPC15', 'MSPTDPC5'}; ...
    };
up.options.opt_config = {'peaks and onsets'; ...
    'nested loops'; ...
    'HRmin=30'; ...
    '20'; ...
    '6'; ...
    %'MSPTDfastv1.1'; ...
    %'Reduced set';...
    };

% % only look at final algorithm
% fields = fieldnames(up.options);
% for field_no = 1 : length(fields)
%     eval(['up.options.' fields{field_no} ' = up.options.' fields{field_no} '(end);'])
% end

end

function generate_results_plots(res, curr_dataset, up)

% setup plots
ftsize = 24;
lwidth = 2;
mksize = 8;
plot_size = [20,20,600,500];
do_quartiles = false;

% go through each potential improvement
for imp_no = 1 : length(up.options.name)

    max_runtime = 0;

    % extract results for this option
    res_rows = [];
    no_options = length(up.options.alg_abbr{imp_no});
    for option_no = 1 : no_options
        res_rows = [res_rows, find(strcmp(res.noQual.num.strategy, up.options.alg_abbr{imp_no}{option_no}))];
    end
    rel_res.f1 = [res.noQual.num.f1_score_med(res_rows), res.noQual.num.f1_score_lq(res_rows), res.noQual.num.f1_score_uq(res_rows)];
    rel_res.runtime = [res.noQual.num.perc_time_med(res_rows), res.noQual.num.perc_time_lq(res_rows), res.noQual.num.perc_time_uq(res_rows)];
    max_runtime = max([max_runtime, max(rel_res.runtime(:,1))]);

    % identify optimal config
    if ~isempty(up.options.opt_config{imp_no})
        opt_config = find(strcmp(up.options.values{imp_no}, up.options.opt_config{imp_no}));
    else
        opt_config = [];
    end
    
    % make plot
    figure('Position', plot_size)

    % - runtime axis
    yyaxis left
    plot(1:no_options, rel_res.runtime(:,1), 'o-b', 'LineWidth', lwidth, 'MarkerSize', mksize - 3)
    hold on
    if do_quartiles
        for q_no = 1:2
            plot(1:no_options, rel_res.runtime(:,1+q_no), '--b', 'LineWidth', lwidth-0.5)
        end
    end
    plot(opt_config, rel_res.runtime(opt_config, 1), 'sb', 'LineWidth', lwidth, 'MarkerSize', mksize+4)
    ylabel('Execution time (%)', 'FontSize', ftsize+2)
    if max_runtime<0.1
        ylim([0.00, ceil(max_runtime*100/2)*2/100])
    else
        ylim([0.00, ceil(max_runtime*100/5)*5/100])
    end
    
    % - f1 score axis
    yyaxis right
    plot(1:no_options, rel_res.f1(:,1), 'o-r', 'LineWidth', lwidth, 'MarkerSize', mksize - 3)
    hold on
    if do_quartiles
        for q_no = 1:2
            plot(1:no_options, rel_res.f1(:,1+q_no), '--r', 'LineWidth', lwidth-0.5)
        end
    end
    plot(opt_config, rel_res.f1(opt_config, 1), 'sr', 'LineWidth', lwidth, 'MarkerSize', mksize+4)
    ylabel('F1-score (%)', 'FontSize', ftsize+2)

    min_yaxis_range = 5;
    max_ylim = max(rel_res.f1(:,1));
    min_ylim = min(rel_res.f1(:,1));
    if (max_ylim-min_ylim)< min_yaxis_range
        med_ylim = mean([max_ylim, min_ylim]);
        max_ylim = med_ylim+(0.5*min_yaxis_range);
        if floor(max_ylim) > max(rel_res.f1(:,1))
            max_ylim = floor(max_ylim);
        elseif ceil(max_ylim-min_yaxis_range) < min(rel_res.f1(:,1))
            max_ylim = ceil(max_ylim);
        end
        if max_ylim > 100 || min_ylim>=(100-min_yaxis_range)
            max_ylim = 100;
        end
        min_ylim = max_ylim-min_yaxis_range;
        ylim([min_ylim, max_ylim])
    end

    xlim([0.5, no_options+0.5])
    set(gca, 'FontSize', ftsize, 'XTick', 1:no_options, 'XTickLabel', up.options.values{imp_no})
    xlabel(up.options.label{imp_no}, 'FontSize', ftsize+2)
    ax = gca;
    ax.YAxis(1).Color = 'b';
    ax.YAxis(2).Color = 'r';
    box off

    rel_dataset_no = find(strcmp(up.datasets, curr_dataset));
    title_txt = up.dataset_names{rel_dataset_no};
    title(title_txt, 'FontSize', ftsize)

    % save plot
    print([up.paths.plot_folder, curr_dataset, '_', up.options.name{imp_no}], '-depsc')
    close all

    % output any text results for this improvement
    if strcmp(up.options.name{imp_no}, 'lmss')
        rel_stat = 'runtime'; rel_or_abs = 'rel';
        orig_config_name = 'peaks and onsets';
        curr_config_name = 'peaks';
        pk_reduction_runtime = calc_reduction(rel_res, rel_stat, imp_no, orig_config_name, curr_config_name, up, rel_or_abs);
        rel_stat = 'f1'; rel_or_abs = 'abs';
        pk_reduction_f1 = calc_reduction(rel_res, rel_stat, imp_no, orig_config_name, curr_config_name, up, rel_or_abs);
        rel_stat = 'runtime'; rel_or_abs = 'rel';
        orig_config_name = 'peaks and onsets';
        curr_config_name = 'onsets';
        onset_reduction_runtime = calc_reduction(rel_res, rel_stat, imp_no, orig_config_name, curr_config_name, up, rel_or_abs);
        rel_stat = 'f1'; rel_or_abs = 'abs';
        onset_reduction_f1 = calc_reduction(rel_res, rel_stat, imp_no, orig_config_name, curr_config_name, up, rel_or_abs);
        fprintf('\nCalculating only a single LMS corresponding to either peaks or onsets resulted in substantial reductions in execution time of %.1f\\%% or %.1f\\%% respectively, compared to calculating LMSs for both peaks and onsets', pk_reduction_runtime, onset_reduction_runtime)
        fprintf('\nThis was accompanied by reductions in \\fscore of %.1f%% or %.1f%%', pk_reduction_f1, onset_reduction_f1);
    elseif strcmp(up.options.name{imp_no}, 'lms_scales')
        rel_stat = 'runtime'; rel_or_abs = 'rel';
        orig_config_name = 'N/2';
        curr_config_name = 'HRmin=40';
        forty_reduction_runtime = calc_reduction(rel_res, rel_stat, imp_no, orig_config_name, curr_config_name, up, rel_or_abs);
        rel_stat = 'f1'; rel_or_abs = 'abs';
        forty_reduction_f1 = calc_reduction(rel_res, rel_stat, imp_no, orig_config_name, curr_config_name, up, rel_or_abs);
        rel_stat = 'runtime'; rel_or_abs = 'rel';
        orig_config_name = 'N/2';
        curr_config_name = 'HRmin=30';
        thirty_reduction_runtime = calc_reduction(rel_res, rel_stat, imp_no, orig_config_name, curr_config_name, up, rel_or_abs);
        rel_stat = 'f1'; rel_or_abs = 'abs';
        thirty_reduction_f1 = calc_reduction(rel_res, rel_stat, imp_no, orig_config_name, curr_config_name, up, rel_or_abs);
        fprintf('\nReducing the number of LMS scales substantially reduced the execution time by %.1f\\%% (with $HR_{min}=40$) or %.1f\\%% (with $HR_{min}=30$).', forty_reduction_runtime, thirty_reduction_runtime)
        fprintf('\nThis was accompanied by reductions in \\fscore of %.1f%% or %.1f%%', forty_reduction_f1, thirty_reduction_f1);
    elseif strcmp(up.options.name{imp_no}, 'fs')
        orig_config_name = 'orig';
        rel_stat = 'runtime'; rel_or_abs = 'rel';
        curr_config_name = '30';
        thirty_reduction_runtime = calc_reduction(rel_res, rel_stat, imp_no, orig_config_name, curr_config_name, up, rel_or_abs);
        rel_stat = 'f1'; rel_or_abs = 'abs';
        thirty_reduction_f1 = calc_reduction(rel_res, rel_stat, imp_no, orig_config_name, curr_config_name, up, rel_or_abs);
        rel_stat = 'runtime'; rel_or_abs = 'rel';
        curr_config_name = '20';
        twenty_reduction_runtime = calc_reduction(rel_res, rel_stat, imp_no, orig_config_name, curr_config_name, up, rel_or_abs);
        rel_stat = 'f1'; rel_or_abs = 'abs';
        twenty_reduction_f1 = calc_reduction(rel_res, rel_stat, imp_no, orig_config_name, curr_config_name, up, rel_or_abs);
        rel_stat = 'runtime'; rel_or_abs = 'rel';
        curr_config_name = '10';
        ten_reduction_runtime = calc_reduction(rel_res, rel_stat, imp_no, orig_config_name, curr_config_name, up, rel_or_abs);
        rel_stat = 'f1'; rel_or_abs = 'abs';
        ten_reduction_f1 = calc_reduction(rel_res, rel_stat, imp_no, orig_config_name, curr_config_name, up, rel_or_abs);
        fprintf('\nReducing the sampling frequency substantially reduced execution time, with values of 30, 20, and 10 Hz reducing execution time by by %.1f\\%%, %.1f\\%% and %.1f\\%% respectively.', thirty_reduction_runtime, twenty_reduction_runtime, ten_reduction_runtime)
        fprintf('\nThese were accompanied by reductions in \\fscore of %.1f\\%%, %.1f\\%% and %.1f\\%%', thirty_reduction_f1, twenty_reduction_f1, ten_reduction_f1);
    elseif strcmp(up.options.name{imp_no}, 'algs')
        orig_config_name = 'MSPTD';
        rel_stat = 'runtime'; rel_or_abs = 'rel';
        curr_config_name = 'MSPTDfastv1.1';
        msptdfast_reduction_runtime = calc_reduction(rel_res, rel_stat, imp_no, orig_config_name, curr_config_name, up, rel_or_abs);
        rel_stat = 'f1'; rel_or_abs = 'abs';
        msptdfast_reduction_f1 = calc_reduction(rel_res, rel_stat, imp_no, orig_config_name, curr_config_name, up, rel_or_abs);
        fprintf('\nThe new \\msptdfast algorithm had a execution time of approximately three-tenths (%.1f\\%%) of the previous \\msptd algorithm', 100-msptdfast_reduction_runtime)
        f1_msptd = rel_res.f1(find(strcmp(up.options.values{imp_no}, orig_config_name)),1);
        f1_msptdfast = rel_res.f1(find(strcmp(up.options.values{imp_no}, curr_config_name)),1);
        fprintf('\nThis was achieved at the expense of a very small reduction in \\fscore of %.1f\\%% (%.1f\\%% vs. %.1f\\%%)', msptdfast_reduction_f1, f1_msptd, f1_msptdfast);
        orig_config_name = 'qppgfast';
        rel_stat = 'runtime'; rel_or_abs = 'rel';
        curr_config_name = 'MSPTDfastv1.1';
        msptdfast_reduction_runtime = calc_reduction(rel_res, rel_stat, imp_no, orig_config_name, curr_config_name, up, rel_or_abs);
        f1_qppgfast = rel_res.f1(find(strcmp(up.options.values{imp_no}, 'qppgfast')),1);
        fprintf('\nIn comparison to \\qppgfast, \\msptdfast had a slightly longer execution time (%.1f\\%% longer) and a comparable \\fscore (%.1f\\%% vs. %.1f\\%%).', -1*msptdfast_reduction_runtime, f1_qppgfast, f1_msptdfast);
        orig_config_name = 'qppgfast';
        rel_stat = 'runtime'; rel_or_abs = 'rel';
        curr_config_name = 'MSPTD';
        msptd_reduction_runtime = calc_reduction(rel_res, rel_stat, imp_no, orig_config_name, curr_config_name, up, rel_or_abs);
        fprintf('\nThe \\msptd algorithm had a %.1f\\%% longer execution time.', -1*msptd_reduction_runtime);

    end

end




end

function red = calc_reduction(rel_res, rel_stat, imp_no, orig_config_name, curr_config_name, up, rel_or_abs)

% identify stat to be used
eval(['stat = rel_res.' rel_stat ';'])

% extract results for the original config and comparator config
orig_config = find(strcmp(up.options.values{imp_no}, orig_config_name));
curr_config = find(strcmp(up.options.values{imp_no}, curr_config_name));

% calculate reduction
if strcmp(rel_or_abs, 'rel')
    red = 100*(stat(orig_config) - stat(curr_config))/stat(orig_config);
elseif strcmp(rel_or_abs, 'abs')
    red = stat(orig_config) - stat(curr_config);
end

end