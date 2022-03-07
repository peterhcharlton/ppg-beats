function assess_multiple_datasets
% ASSESS_MULTIPLE_DATASETS  Assess on multiple datasets.
%   ASSESS_MULTIPLE_DATASETS assesses the performance of multiple
%   PPG beat detectors across multiple datasets.
%   
%   # Inputs
%   * none
%   * NB: the 'ppg_beat_detector_assessment' function, which is called by this function, requires a dataset to analyse
%   
%   # Outputs
%   * Plots illustrating the performance of PPG beat detectors.
%   * NB: the 'ppg_beat_detector_assessment' function, which is called by this function, creates its own outputs too
%   
%   # Exemplary usage
%   
%       run_ppg_beat_detector_assessment      % runs the PPG beat detector assessment
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

%% Setup
% This contains settings (including filepaths) which must be set before running the script.
up = setup_universal_params;

% Specify options for beat detector assessment
options = specify_options;

%% Perform analysis on each dataset in turn

do_analysis = 1;
if do_analysis
    for dataset_no = 1 : length(up.datasets)
        
        curr_dataset = up.datasets{dataset_no};
        
        % Specify the path of the Matlab data file for this dataset
        options.dataset_file = specify_path_of_dataset_file(curr_dataset);

        % Assess performance of PPG beat detectors on this dataset
        assess_beat_detectors(curr_dataset, options);
        
    end
end

%% Retrieve results for strategies, and provide latex table of results of selected strategy across all datasets
selected_strategy = 'MSPTD__accel_comb_MSPTD_qppg';
selected_strategy = 'MSPTD__none';
vars = {'ppv', 'sens', 'f1_score', 'acc_ibi', 'prop_ibi', 'acc_hr', 'prop_hr'};
types = {'med', 'lq', 'uq'};
fprintf('\n - Performance of selected strategy: ')
[total_pts, total_durn] = deal(0);
for dataset_no = 1 : length(up.datasets)
    
    curr_dataset = up.datasets{dataset_no};
    fprintf('\n%s', create_title_text(curr_dataset))
    
    % Specify options for this dataset
    options = specify_options;
    options.dataset_file = specify_path_of_dataset_file(curr_dataset);
    
    % Specify processing folder for this dataset
    [folder_containing_raw_data_file, ~, ~] = fileparts(options.dataset_file);
    eval(sprintf('up.paths.processing_folder = ''%s%s%s%s%s'';', folder_containing_raw_data_file, filesep, 'proc_data_', curr_dataset, filesep));
    
    % Load all results for this dataset
    res_file = [up.paths.processing_folder, 'ppg_detect_perf.mat'];
    load(res_file, 'ppg_strategy_perf'); clear res_file
    
    % store total beats and total duration for this dataset
    if ~strcmp(curr_dataset, 'mimic_test_n') && ~strcmp(curr_dataset, 'mimic_test_a')
        eval(['temp = ppg_strategy_perf.raw.' selected_strategy '.durn_total;']);
        total_pts = total_pts + length(temp);
        total_durn = total_durn + sum(temp(~isnan(temp)));
    end
    
    % Extract raw stats for each relevant strategy
    rel_strategy_els = find(contains(ppg_strategy_perf.stats.strategies, '__none'));
    raw_res.rel_strategies = ppg_strategy_perf.stats.strategies(rel_strategy_els);
    vars = {'f1_score', 'sens', 'ppv', 'durn_both', 'no_beats'};
    for curr_var_no = 1 : length(vars)
        curr_var = vars{curr_var_no};
        for strategy_no = 1 : length(rel_strategy_els)
            curr_strategy = ppg_strategy_perf.stats.strategies{rel_strategy_els(strategy_no)};
            eval(['temp = ppg_strategy_perf.raw.' curr_strategy '.' curr_var ';']);
            eval(['raw_res.' curr_dataset '.' curr_var '(strategy_no,:) = temp;']);
        end
    end
    clear temp curr_strategy strategy_no curr_var rel_strategy_els
    
    % Extract results for the selected strategy
    strategy_el = find(strcmp(ppg_strategy_perf.stats.strategies, selected_strategy));
    if isempty(strategy_el)
        strategy_el = find(strcmp(ppg_strategy_perf.stats.strategies, strrep(selected_strategy, '_accel', '')));
    end
    if ~isempty(strategy_el)
       for var_no = 1 : length(vars)
            curr_var = vars{var_no};
            
            % - variables requiring median and quartile ranges
            if ~sum(strcmp(curr_var, {'no_beats'}))
                for type_no = 1 : length(types)
                    curr_type = types{type_no};
                    eval(['curr_temp = ppg_strategy_perf.stats.' curr_var '.' curr_type '(strategy_el);']);
                    if strcmp(curr_var, 'durn_both')
                        curr_temp = curr_temp./60;
                    end
                    eval(['temp.' curr_type ' = curr_temp;']);
                end
                fprintf(' & %.1f (%.1f - %.1f)', temp.med, temp.lq, temp.uq)
                
                % - variables requiring totalling
            else
                curr_strategy = ppg_strategy_perf.stats.strategies{strategy_el};
                eval(['curr_temp_data = ppg_strategy_perf.raw.' curr_strategy '.' curr_var ';']);
                temp.total = sum(curr_temp_data(~isnan(curr_temp_data)));
                fprintf(' & %d', temp.total)                
            end
            
            % - store results for this variable on this dataset
            eval(['dataset_results.' curr_var ' = temp;']); clear temp
        end
        fprintf(' \\\\')
        
        % store results
        eval(['strategy_res.' curr_dataset ' = dataset_results;']);
    end
    
    clear dataset_results ppg_strategy_perf options strategy_el var_no type_no curr_var curr_type
    
end
fprintf('\n\n - Total of %d subjects (although there is some overlap between datasets) and %.1f hrs of continuous recording\n\n', total_pts, total_durn/(60*60))
clear vars types

%% Collate results tables for all datasets
for dataset_no = 1 : length(up.datasets)
    
    curr_dataset = up.datasets{dataset_no};
    
    % Specify options for this dataset
    options = specify_options;
    options.dataset_file = specify_path_of_dataset_file(curr_dataset);
    
    % Specify processing folder for this dataset
    [folder_containing_raw_data_file, ~, ~] = fileparts(options.dataset_file);
    eval(sprintf('up.paths.processing_folder = ''%s%s%s%s%s'';', folder_containing_raw_data_file, filesep, 'proc_data_', up.datasets{dataset_no}, filesep));
    
    % Load summary results for this dataset
    res_file = [up.paths.processing_folder, 'ppg_detect_stats.mat'];
    dataset_res = load(res_file, 'res'); clear res_file
    
    % Create results tables
    up.analysis.beat_detectors = options.beat_detectors;
    dataset_res = dataset_res.res;
    
    % store results
    eval(['res.' up.datasets{dataset_no} ' = dataset_res;']);
    clear dataset_res ppg_strategy_perf options
    
end

%% Generate comparison figures

do_comparison_figs = 1;
if do_comparison_figs
    up.paths.plots_folder = [up.paths.plots_root_folder, 'other', filesep];
    create_comparison_figures(up.datasets, raw_res, res, up);
end

%% Generate results figures

fprintf('\n\n ~~~ Making performance figures ~~~')
for do_training = [1,0]
    if do_training == 1
        fprintf('\n - Training datasets')
        curr_datasets = up.training_datasets;
        up.paths.plots_folder = [up.paths.plots_root_folder, 'training', filesep];  
    elseif do_training == 0
        fprintf('\n - Testing datasets')
        curr_datasets = up.testing_datasets; 
        up.paths.plots_folder = [up.paths.plots_root_folder, 'validation', filesep];  
    else
        curr_datasets = {'ppg_dalia_sitting', 'ppg_dalia_working', 'ppg_dalia_walking'};
        up.paths.plots_folder = [up.paths.plots_root_folder, 'other', filesep];  
    end
    
    if ~exist(up.paths.plots_folder, 'dir'), mkdir(up.paths.plots_folder), end
        
    create_perf_box_plots(res, up, curr_datasets, do_training);
end

end

function up = setup_universal_params

% - Datasets to be analysed
up.training_datasets = {'capnobase', 'bidmc', 'mimic_train_all', 'mimic_test_all', 'wesad_meditation', 'wesad_amusement', 'wesad_baseline', 'wesad_stress', 'ppg_dalia_sitting', 'ppg_dalia_working', 'ppg_dalia_cycling', 'ppg_dalia_walking'}; 
% up.training_datasets = {'mimic_non_af', 'mimic_af', 'ppg_dalia_sitting', 'ppg_dalia_working', 'ppg_dalia_walking'}; 
up.testing_datasets = {'mimic_B', 'mimic_W', 'mimic_test_a', 'mimic_test_n', 'mimic_non_af', 'mimic_af'};
up.datasets = [up.testing_datasets, up.training_datasets];

% Specify path at which to save plots
up.paths.plots_root_folder = '/Users/petercharlton/Google Drive/Work/Publications/In Preparation/2021 Assess PPG beat detectors/figures/';

% Specify comparisons
up.settings.comparisons = {'mimic_B', 'mimic_W'; 'mimic_test_a', 'mimic_test_n'; 'mimic_non_af', 'mimic_af'};

%% Display startup message
display_startup_message;

end

function display_startup_message
% Displays a startup message, including details of the licence

licence_details = ['\n\n RUN_PPG_BEAT_DETECTOR_ASSESSMENT  Copyright (C) 2022  Peter H. Charlton',...
    '\n This program comes with ABSOLUTELY NO WARRANTY', ... 
    '\n and is available under the GNU public license.', ...
    '\n For details see the accompanying LICENSE.txt file.\n'];

fprintf(licence_details)

end

function create_comparison_figures(datasets, raw_res, res, uParams)

fprintf('\n ~~~ Making comparison figures ~~~')

% setup
ftsize = 22;

if sum(strcmp(datasets, 'ppg_dalia_working'))
    primary_dataset = 'ppg_dalia_working';
elseif sum(strcmp(datasets, 'wesad_baseline'))
    primary_dataset = 'wesad_baseline';
else
    primary_dataset = datasets{1};
end
eval(['orig_beat_detector_names = res.' primary_dataset '.noQual.num.strategy;']);
orig_beat_detector_names = orig_beat_detector_names(~strcmp(orig_beat_detector_names, 'qppg'));
beat_detector_names = correct_beat_detector_names(orig_beat_detector_names);

% identify relevant results set
curr_strategy_set = 'noQual';

% cycle through each comparison (i.e. pair of datasets)
no_comparisons = size(uParams.settings.comparisons,1);
for comparison_no = 1 : no_comparisons
    
    fprintf('\n - Making comparison figure for %s and %s: ', uParams.settings.comparisons{comparison_no,1}, uParams.settings.comparisons{comparison_no,2})
    
    % cycle through each statistic
    eval(['vars = fieldnames(raw_res.' uParams.settings.comparisons{comparison_no,1} ');']);
    vars = vars(~strcmp(vars, 'rel_strategies') & ~strcmp(vars, 'durn_ibi') & ~strcmp(vars, 'durn_both') & ~strcmp(vars, 'no_beats'));
    for var_no = 1 : length(vars)
        curr_var = vars{var_no};
        
        fprintf('%s, ', curr_var)
        
        figure('Position', [20,20,600,450]);
        subplot('Position', [0.1,0.1,0.89,0.8])

        % extract subject-level results for each dataset
        eval(['res1.raw.v = raw_res.' uParams.settings.comparisons{comparison_no,1} '.' curr_var ';']);
        res1.raw.strategies = raw_res.rel_strategies;
        eval(['res2.raw.v = raw_res.' uParams.settings.comparisons{comparison_no,2} '.' curr_var ';']);
        res2.raw.strategies = raw_res.rel_strategies;
        
        strategies_log = ~contains(raw_res.rel_strategies, 'qppg_');
        
        % setup boxplot
        group_no = 0; [curr_res, groups, positions, box_nos, pc10s, pc90s, xticks] = deal([]);
        
        % cycle through each beat detection strategy
        for beat_detector_no = 1 : length(beat_detector_names)
            curr_beat_detector = orig_beat_detector_names{beat_detector_no};
            res1_raw_strategy_no = find(strcmp(res1.raw.strategies, [curr_beat_detector, '__none']));
            res2_raw_strategy_no = find(strcmp(res2.raw.strategies, [curr_beat_detector, '__none']));
            
            % Test whether the use of a particular beat detection strategy on these two datasets resulted in significantly different results.
            temp_p(comparison_no,beat_detector_no,var_no) = ranksum(res1.raw.v(res1_raw_strategy_no,:), res2.raw.v(res2_raw_strategy_no,:));
            
            % Store results to make boxplot
            curr_res = [curr_res, res1.raw.v(res1_raw_strategy_no,:), res2.raw.v(res2_raw_strategy_no,:)];
            groups = [groups, group_no*ones(1,size(res1.raw.v,2)), (group_no+1)*ones(1,size(res2.raw.v,2))];
            positions = [positions, group_no, (group_no+0.75)];
            box_nos = [box_nos, (group_no+1), (group_no+2)];
            xticks = [xticks, (group_no+(0.75/2)), (group_no+(0.75/2))];
            pc10s = [pc10s, quantile(res1.raw.v(res1_raw_strategy_no,:),0.1), quantile(res2.raw.v(res2_raw_strategy_no,:),0.1)];
            pc90s = [pc90s, quantile(res1.raw.v(res1_raw_strategy_no,:),0.9), quantile(res2.raw.v(res2_raw_strategy_no,:),0.9)];
            group_no = group_no+2;
            
        end
        clear strategy_no
        
        % create box plot
        boxplot(curr_res, groups, 'positions', positions, 'MedianStyle', 'target', 'Symbol', 'r', 'Widths', 0.6);
        set(gca, 'XTickLabelRotation', 90, 'FontSize', ftsize-2)
        set(gca,'XGrid', 'on', 'YGrid', 'on')
        box off
        beat_detector_xticks = unique(xticks);
        set(gca, 'XTick', beat_detector_xticks, 'XTickLabel', beat_detector_names') % used to be beat_detectors
        [ylims, ylab_txt] = get_ylabel_settings(curr_var);
        ylim([65,100]);
        ylabel(ylab_txt, 'FontSize', ftsize)
        
        % adjust whiskers to be 10th and 90th percentiles, inspired by: https://www.mathworks.com/matlabcentral/fileexchange/22526-box-plot-with-whiskers-plotted-at-fixed-percentiles
        lower_whiskers = findobj(gca,'Tag','Lower Whisker');
        lower_whisker_ends = findobj(gca,'Tag','Lower Adjacent Value');
        upper_whiskers = findobj(gca,'Tag','Upper Whisker');
        upper_whisker_ends = findobj(gca,'Tag','Upper Adjacent Value');
        for whisker_no = 1 : length(lower_whiskers)
            rel_el = find(positions == lower_whiskers(whisker_no).XData(1),1);
            set(lower_whiskers(whisker_no),'YData',[pc10s(rel_el), lower_whiskers(whisker_no).YData(2)], 'LineStyle', '-');
            set(lower_whisker_ends(whisker_no),'YData',ones(1,2)*pc10s(rel_el));
            set(upper_whiskers(whisker_no),'YData',[upper_whiskers(whisker_no).YData(1), pc90s(rel_el)], 'LineStyle', '-');
            set(upper_whisker_ends(whisker_no),'YData',ones(1,2)*pc90s(rel_el));
        end
        
        % add colour to boxes, based on: https://www.mathworks.com/matlabcentral/answers/392679-how-to-fill-boxes-in-boxplot-with-different-colors
        h = findobj(gca,'Tag','Box');
        for j=1:length(h)
            % if the first box in a pair
            if rem(ceil(median(h(j).XData)),2) == 0
                face_color = 0.4*ones(1,3); % greyscale
                face_color = [1,0,0];
                hleg(1,1) = patch(get(h(j),'XData'),get(h(j),'YData'),face_color,'FaceAlpha',.5);
            else
                face_color = 0.8*ones(1,3); % greyscale
                face_color = [0,0,1];
                hleg(1,2) = patch(get(h(j),'XData'),get(h(j),'YData'),face_color,'FaceAlpha',.5);
            end
        end
        
        % change color of median target point
        h1 = findobj(gca,'Tag','MedianInner');
        h2 = findobj(gca,'Tag','MedianInner');
        for j=1:length(h1)
            hold on
            plot(h1(j).XData, h1(j).YData, 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k')
            plot(h2(j).XData, h2(j).YData, 'ok', 'MarkerSize', 10)
            delete(h1(j))
            delete(h2(j))
        end
        
        % create legend
        if isequal(uParams.settings.comparisons(comparison_no,:), {'mimic_B','mimic_W'})
            leg_labels = {'Black', 'White'};
        elseif isequal(uParams.settings.comparisons(comparison_no,:), {'mimic_test_a','mimic_test_n'})
            leg_labels = {'Adults', 'Neonates'};
        elseif isequal(uParams.settings.comparisons(comparison_no,:), {'mimic_non_af','mimic_af'})
            leg_labels = {'Non-AF', 'AF'};
        end
        legend(hleg, leg_labels, 'Position', [0.33,0.9,0.5,0.1], 'Orientation', 'Horizontal')
        
        % add significant differences
        pos = get(gca, 'Position');
        xlims = get(gca, 'XLim');
        holm_sidak_alpha = 0.0014;
        no_sig = 0;
        for strategy_no = 1 :size(temp_p,2)
            if strcmp(vars{var_no}, 'f1_score') && temp_p(comparison_no,strategy_no,var_no) < holm_sidak_alpha
                pos_x_coord = pos(1) + pos(3)*((beat_detector_xticks(strategy_no)-xlims(1))/(xlims(2)-xlims(1)));
                y_coord = .28-0.05*rem(no_sig,2);
                dim = [pos_x_coord-0.05 y_coord .1 .1];
                str = 'p<0.002';
                annotation('textbox',dim,'String',str,'FitBoxToText','on','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize', ftsize-8, 'LineStyle', 'none');
                no_sig = no_sig + 1;
            elseif temp_p(comparison_no,strategy_no,var_no) < 0.05
                pos_x_coord = pos(1) + pos(3)*((beat_detector_xticks(strategy_no)-xlims(1))/(xlims(2)-xlims(1)));
                y_coord = .28-0.05*rem(no_sig,2);
                dim = [pos_x_coord-0.05 y_coord .1 .1];
                str = 'p<0.05';
                annotation('textbox',dim,'String',str,'FitBoxToText','on','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize', ftsize-8, 'LineStyle', 'none');
                no_sig = no_sig + 1;
            end
            
        end
        clear res1 res2
        
        % make slightly wider to include any p-value annotations
        xlim([xlims(1), xlims(2)+1])
        
        % add x-axis label
        xlabel('Beat detector', 'FontSize', ftsize)
        
        % save figure
        fig_name = ['comparison_' leg_labels{1} '_' leg_labels{2} '_' curr_var];
        save_figure(gcf, fig_name, uParams);
        
        close all;
        
    end
    clear var_no curr_var
    
end
clear comparison_no

% find significant differences
alpha = 0.05;
temp_p_f1_score = temp_p(:,:,strcmp(vars,'f1_score'));
all_p = sort(temp_p_f1_score(:));
temp_sig_diffs = zeros(length(all_p),1);
for comparison_no = 1 : length(all_p)
    % Holm's correction for only the number of remaining comparisons
    no_tests = length(all_p)-comparison_no+1;
    % Sidak's correction for the number of remaining comparisons
    alpha_sidak = 1 - ((1-alpha)^(1/no_tests));
    if all_p(comparison_no) < alpha_sidak
        temp_sig_diffs(comparison_no) = 1;
    else
        break
    end
end
curr_thresh = round(10000*alpha_sidak)/10000;
if ~isequal(curr_thresh, holm_sidak_alpha)
    warning('Alpha value needs updating')
end

end

function beat_detectors = correct_beat_detector_names(beat_detectors)
beat_detectors = strrep(beat_detectors, 'qppgfast', 'qppg');
beat_detectors = strrep(beat_detectors, 'PPGPulses', 'Pulses');
end

function res = create_results_tables(ppg_strategy_perf, uParams)

%% Create rankings

% - identify performance metrics
eval(['perf_metrics = fieldnames(ppg_strategy_perf.raw.' ppg_strategy_perf.stats.strategies{1} ');']);

% --- Best strategies with no quality assessment
rel_strategies = contains(ppg_strategy_perf.stats.strategies, '__none');
[~,order] = sort(ppg_strategy_perf.stats.f1_score.med(rel_strategies), 'descend');
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
        temp = strfind(strategy{el_no}, '__');
        strategy{el_no} = strategy{el_no}{1}(1:temp{1}-1); clear temp
        clear med lq uq curr_txt
    end
    clear el_no curr_metric
end
clear metric_no els
res.no_qual = table(strategy, sens, ppv, f1_score, new_score);
clear strategy sens ppv f1_score new_score

% --- Best strategies with accel only
rel_strategies = false(size(ppg_strategy_perf.stats.strategies));
for s = 1 : length(rel_strategies)
    if strcmp(ppg_strategy_perf.stats.strategies{s}(end-6:end), '__accel')
        rel_strategies(s) = true;
    end
end
% - skip if there are no 'accel' strategies
if sum(rel_strategies) ~= 0
    [~,order] = sort(ppg_strategy_perf.stats.f1_score.med(rel_strategies), 'descend');
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
            temp = strfind(strategy{el_no}, '__');
            strategy{el_no} = strategy{el_no}{1}(1:temp{1}-1); clear temp
            clear med lq uq curr_txt
        end
        clear el_no curr_metric
    end
    clear metric_no els
    res.accel = table(strategy, sens, ppv, f1_score, new_score);
    clear strategy sens ppv f1_score new_score
end

% --- Best strategies without accel with qual assessment using same beat detector
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
    [~,order] = sort(ppg_strategy_perf.stats.f1_score.med(temp_rel_strategies), 'descend');
    rel_strategies(temp_rel_strategies(order(1))) = true;
end
%  - part 2 : sort these strategies
[~,order] = sort(ppg_strategy_perf.stats.f1_score.med(rel_strategies), 'descend');
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

% --- Best strategies with accel with qual assessment using same beat detector
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
    [~,order] = sort(ppg_strategy_perf.stats.f1_score.med(temp_rel_strategies), 'descend');
    if ~isempty(order)
        rel_strategies(temp_rel_strategies(order(1))) = true;
    end
end
% - skip if there are no 'accel' strategies
if sum(rel_strategies) ~= 0
    %  - part 2 : sort these strategies
    [~,order] = sort(ppg_strategy_perf.stats.f1_score.med(rel_strategies), 'descend');
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
    res.comb_accel = table(strategy, sens, ppv, f1_score, new_score);
    clear strategy sens ppv f1_score new_score
end

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
[~,order] = sort(ppg_strategy_perf.stats.f1_score.med(rel_strategies), 'descend');
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

function options = specify_options
% Specify options for analysis

% Setup default options
options = struct;

% Specify the beat detectors to be used
options.beat_detectors = {'MSPTD', 'qppg', 'ABD', 'COppg'};
options.beat_detectors = {'SPAR3', 'SPAR7', 'IMS', 'AMPD', 'MSPTD', 'ABD', 'qppg', 'HeartPy', 'COppg', 'Pulses'};
options.beat_detectors = {'SWT', 'SPAR', 'IMS', 'AMPD', 'MSPTD', 'ABD', 'qppg', 'qppgfast', 'HeartPy', 'COppg', 'PPGPulses', 'ERMA', 'PWD', 'PDA', 'WFD'};
% Extras: 'ATmin', 'ATmax', 

% Specify the downsampling strategy
options.do_downsample = 1;     % downsample PPG signals
options.downsample_freq = 100; % frequency in Hz

close all
end

function dataset_file = specify_path_of_dataset_file(dataset)
% Specify the file path of the Matlab data file containing the data for this dataset.

switch dataset
    case 'mimic_af'
        dataset_file = '/Users/petercharlton/Documents/Data/mimiciii_ppg_af_beat_detection/mimic_af_data.mat';
    case 'mimic_non_af'
        dataset_file = '/Users/petercharlton/Documents/Data/mimiciii_ppg_af_beat_detection/mimic_non_af_data.mat';
    case 'ppg_dalia_sitting'
        dataset_file = '/Users/petercharlton/Documents/Data/ppg_dalia/conv_data/ppg_dalia_sitting_data.mat';
    case 'ppg_dalia_working'
        dataset_file = '/Users/petercharlton/Documents/Data/ppg_dalia/conv_data/ppg_dalia_working_data.mat';
    case 'ppg_dalia_walking'
        dataset_file = '/Users/petercharlton/Documents/Data/ppg_dalia/conv_data/ppg_dalia_walking_data.mat';
    case 'ppg_dalia_lunch_break'
        dataset_file = '/Users/petercharlton/Documents/Data/ppg_dalia/conv_data/ppg_dalia_lunch_break_data.mat';
    case 'ppg_dalia_stair_climbing'
        dataset_file = '/Users/petercharlton/Documents/Data/ppg_dalia/conv_data/ppg_dalia_stair_climbing_data.mat';
    case 'ppg_dalia_car_driving'
        dataset_file = '/Users/petercharlton/Documents/Data/ppg_dalia/conv_data/ppg_dalia_car_driving_data.mat';
    case 'ppg_dalia_table_soccer'
        dataset_file = '/Users/petercharlton/Documents/Data/ppg_dalia/conv_data/ppg_dalia_table_soccer_data.mat';
    case 'ppg_dalia_cycling'
        dataset_file = '/Users/petercharlton/Documents/Data/ppg_dalia/conv_data/ppg_dalia_cycling_data.mat';
    case 'capnobase'
    	dataset_file = '/Users/petercharlton/Documents/Data/capnobase_2021/conv_data/capnobase_data.mat';
	case 'wesad_amusement'
    	dataset_file = '/Users/petercharlton/Documents/Data/WESAD/conv_data/wesad_amusement_data.mat';
	case 'wesad_baseline'
    	dataset_file = '/Users/petercharlton/Documents/Data/WESAD/conv_data/wesad_baseline_data.mat';
	case 'wesad_meditation'
    	dataset_file = '/Users/petercharlton/Documents/Data/WESAD/conv_data/wesad_meditation_data.mat';
	case 'wesad_stress'
    	dataset_file = '/Users/petercharlton/Documents/Data/WESAD/conv_data/wesad_stress_data.mat';
    case 'bidmc'
        dataset_file = '/Users/petercharlton/Documents/Data/BIDMC_dataset/bidmc_data.mat';
    case 'mimic_B'
        dataset_file = '/Users/petercharlton/Documents/Data/mimiciii_ppg_ethnicity_beat_detection/mimic_B_data.mat';
    case 'mimic_W'
        dataset_file = '/Users/petercharlton/Documents/Data/mimiciii_ppg_ethnicity_beat_detection/mimic_W_data.mat';
    case 'mimic_a'
        dataset_file = '/Users/petercharlton/Documents/Data/mimiciii_ppg_neonate_adult_beat_detection/mimic_a_data.mat';
    case 'mimic_n'
        dataset_file = '/Users/petercharlton/Documents/Data/mimiciii_ppg_neonate_adult_beat_detection/mimic_n_data.mat';
    case 'mimic_train_all'
        dataset_file = '/Users/petercharlton/Documents/Data/mimiciii_ppg_neonate_adult_beat_detection/mimic_train_all_data.mat';
    case 'mimic_test_all'
        dataset_file = '/Users/petercharlton/Documents/Data/mimiciii_ppg_neonate_adult_beat_detection/mimic_test_all_data.mat';
    case 'mimic_test_a'
        dataset_file = '/Users/petercharlton/Documents/Data/mimiciii_ppg_neonate_adult_beat_detection/mimic_test_a_data.mat';
    case 'mimic_test_n'
        dataset_file = '/Users/petercharlton/Documents/Data/mimiciii_ppg_neonate_adult_beat_detection/mimic_test_n_data.mat';
end

end

function create_perf_box_plots(res, uParams, curr_datasets, do_training)

if do_training ~= 0
    %% Beat detection box plot
    datasets = fieldnames(res);
    datasets = curr_datasets;
    if sum(strcmp(datasets, 'ppg_dalia_working'))
        primary_dataset = 'ppg_dalia_working';
    elseif sum(strcmp(datasets, 'wesad_baseline'))
        primary_dataset = 'wesad_baseline';
    else
        primary_dataset = datasets{1};
    end
    eval(['beat_detectors = res.' primary_dataset '.noQual.num.strategy;']);
    curr_beat_detectors = beat_detectors(~strcmp(beat_detectors, 'qppg'));
    rel_metric = 'f1_score';
    rel_strategy_set = 'noQual';
    fig_name = 'beat_detector_f1_boxplot';
    curr_datasets = datasets(~strcmp(datasets, 'ppg_dalia_walking'));
    curr_datasets = datasets;
    make_boxplot_figure(res, primary_dataset, curr_datasets, curr_beat_detectors, rel_metric, rel_strategy_set, fig_name, uParams);
    
    %% Beat detection PPV box plot
    datasets = fieldnames(res);
    datasets = curr_datasets;
    if sum(strcmp(datasets, 'ppg_dalia_working'))
        primary_dataset = 'ppg_dalia_working';
    else
        primary_dataset = datasets{1};
    end
    eval(['beat_detectors = res.' primary_dataset '.noQual.num.strategy;']);
    curr_beat_detectors = beat_detectors(~strcmp(beat_detectors, 'qppg'));
    rel_metric = 'ppv';
    rel_strategy_set = 'noQual';
    fig_name = 'beat_detector_ppv_boxplot';
    curr_datasets = datasets(~strcmp(datasets, 'ppg_dalia_walking'));
    curr_datasets = datasets;
    make_boxplot_figure(res, primary_dataset, curr_datasets, curr_beat_detectors, rel_metric, rel_strategy_set, fig_name, uParams);
    
    %% Beat detection sens box plot
    datasets = fieldnames(res);
    datasets = curr_datasets;
    if sum(strcmp(datasets, 'ppg_dalia_working'))
        primary_dataset = 'ppg_dalia_working';
    else
        primary_dataset = datasets{1};
    end
    eval(['beat_detectors = res.' primary_dataset '.noQual.num.strategy;']);
    curr_beat_detectors = beat_detectors(~strcmp(beat_detectors, 'qppg'));
    rel_metric = 'sens';
    rel_strategy_set = 'noQual';
    fig_name = 'beat_detector_sens_boxplot';
    curr_datasets = datasets(~strcmp(datasets, 'ppg_dalia_walking'));
    curr_datasets = datasets;
    make_boxplot_figure(res, primary_dataset, curr_datasets, curr_beat_detectors, rel_metric, rel_strategy_set, fig_name, uParams);
    
    %% HRs no qual box plot
    datasets = fieldnames(res);
    if sum(strcmp(datasets, 'ppg_dalia_working'))
        primary_dataset = 'ppg_dalia_working';
    else
        primary_dataset = datasets{1};
    end
    rel_metric = 'acc_hr';
    rel_strategy_set = 'noQual';
    fig_name = 'HR_no_qual_boxplot';
    make_boxplot_figure(res, primary_dataset, curr_datasets, curr_beat_detectors, rel_metric, rel_strategy_set, fig_name, uParams);
    rel_metric = 'mae_hr';
    rel_strategy_set = 'noQual';
    fig_name = 'HR_mae_no_qual_boxplot';
    make_boxplot_figure(res, primary_dataset, curr_datasets, curr_beat_detectors, rel_metric, rel_strategy_set, fig_name, uParams);
    
    %% IBIs no qual box plot
    datasets = fieldnames(res);
    if sum(strcmp(datasets, 'ppg_dalia_working'))
        primary_dataset = 'ppg_dalia_working';
    else
        primary_dataset = datasets{1};
    end
    rel_metric = 'acc_ibi';
    rel_strategy_set = 'noQual';
    fig_name = 'IBI_acc_no_qual_boxplot';
    make_boxplot_figure(res, primary_dataset, curr_datasets, curr_beat_detectors, rel_metric, rel_strategy_set, fig_name, uParams);
    rel_metric = 'mae_ibi';
    rel_strategy_set = 'noQual';
    fig_name = 'IBI_mae_no_qual_boxplot';
    make_boxplot_figure(res, primary_dataset, curr_datasets, curr_beat_detectors, rel_metric, rel_strategy_set, fig_name, uParams);
    
end

if do_training ~= 1
    
    %% Beat detection with and without accel only and in combination qual box plot
    datasets = fieldnames(res);
    datasets = curr_datasets;
    if sum(strcmp(datasets, 'wesad_baseline'))
        primary_dataset = 'wesad_baseline';
    elseif sum(strcmp(datasets, 'ppg_dalia_working'))
        primary_dataset = 'ppg_dalia_working';
    else
        primary_dataset = datasets{1};
    end
    eval(['beat_detectors = res.' primary_dataset '.noQual.num.strategy;']);
    curr_beat_detectors = beat_detectors(~strcmp(beat_detectors, 'SPAR3'));
    rel_metric = 'ppv';
    %curr_datasets = datasets(~cellfun(@isempty, strfind(datasets, 'ppg_dalia')));
    %curr_datasets = datasets(~cellfun(@isempty, strfind(datasets, 'wesad')));
    rel_strategy_sets = {'noQual', 'accel', 'comb_noAccel', 'comb_Accel'};
    fig_name = 'beat_detector_ppv_with_without_accel_boxplot';
    make_boxplot_figure_multiple_approaches(res, primary_dataset, curr_datasets, curr_beat_detectors, rel_metric, rel_strategy_sets, fig_name, uParams);
    rel_metric = 'sens';
    %curr_datasets = datasets(~cellfun(@isempty, strfind(datasets, 'ppg_dalia')));
    %curr_datasets = datasets(~cellfun(@isempty, strfind(datasets, 'wesad')));
    fig_name = 'beat_detector_sens_with_without_accel_boxplot';
    make_boxplot_figure_multiple_approaches(res, primary_dataset, curr_datasets, curr_beat_detectors, rel_metric, rel_strategy_sets, fig_name, uParams);
    rel_metric = 'f1_score';
    %curr_datasets = datasets(~cellfun(@isempty, strfind(datasets, 'ppg_dalia')));
    %curr_datasets = datasets(~cellfun(@isempty, strfind(datasets, 'wesad')));
    fig_name = 'beat_detector_f1_with_without_accel_boxplot';
    make_boxplot_figure_multiple_approaches(res, primary_dataset, curr_datasets, curr_beat_detectors, rel_metric, rel_strategy_sets, fig_name, uParams);
    
    %% HRs with and without accel only and in combination qual box plot
    datasets = fieldnames(res);
    if sum(strcmp(datasets, 'ppg_dalia_working'))
        primary_dataset = 'ppg_dalia_working';
    else
        primary_dataset = datasets{1};
    end
    rel_metric = 'mae_hr';
    %curr_datasets = datasets(~cellfun(@isempty, strfind(datasets, 'ppg_dalia')));
    %curr_datasets = datasets(~cellfun(@isempty, strfind(datasets, 'wesad')));
    fig_name = 'HR_mae_with_without_accel_boxplot';
    make_boxplot_figure_multiple_approaches(res, primary_dataset, curr_datasets, curr_beat_detectors, rel_metric, rel_strategy_sets, fig_name, uParams);
    
    %% IBIs with and without accel only and in combination qual box plot
    datasets = fieldnames(res);
    if sum(strcmp(datasets, 'ppg_dalia_working'))
        primary_dataset = 'ppg_dalia_working';
    else
        primary_dataset = datasets{1};
    end
    rel_metric = 'mae_ibi';
    %curr_datasets = datasets(~cellfun(@isempty, strfind(datasets, 'ppg_dalia')));
    %curr_datasets = datasets(~cellfun(@isempty, strfind(datasets, 'wesad')));
    fig_name = 'IBI_mae_with_without_accel_boxplot';
    make_boxplot_figure_multiple_approaches(res, primary_dataset, curr_datasets, curr_beat_detectors, rel_metric, rel_strategy_sets, fig_name, uParams);
    
    %% HRs with and without accel only and in combination qual box plot
    datasets = fieldnames(res);
    if sum(strcmp(datasets, 'ppg_dalia_working'))
        primary_dataset = 'ppg_dalia_working';
    else
        primary_dataset = datasets{1};
    end
    rel_metric = 'acc_hr';
    %curr_datasets = datasets(~cellfun(@isempty, strfind(datasets, 'ppg_dalia')));
    %curr_datasets = datasets(~cellfun(@isempty, strfind(datasets, 'wesad')));
    fig_name = 'HR_acc_with_without_accel_boxplot';
    make_boxplot_figure_multiple_approaches(res, primary_dataset, curr_datasets, curr_beat_detectors, rel_metric, rel_strategy_sets, fig_name, uParams);
    
end

return

%% IBIs with and without accel only qual box plot
datasets = fieldnames(res);
if sum(strcmp(datasets, 'ppg_dalia_working'))
    primary_dataset = 'ppg_dalia_working';
else
    primary_dataset = datasets{1};
end
rel_metric = 'acc_ibi';
rel_strategy_sets = {'noQual', 'accel'};
fig_name = 'IBI_acc_with_without_accel_only_boxplot';
curr_datasets = datasets(~cellfun(@isempty, strfind(datasets, 'ppg_dalia')));
%curr_datasets = datasets(~cellfun(@isempty, strfind(datasets, 'wesad')));
make_boxplot_figure_two_approaches(res, primary_dataset, curr_datasets, curr_beat_detectors, rel_metric, rel_strategy_sets, fig_name, uParams);
rel_strategy_sets = {'noQual', 'accel', 'comb_noAccel'};
fig_name = 'IBI_acc_with_without_accel_comb_boxplot';
make_boxplot_figure_multiple_approaches(res, primary_dataset, curr_datasets, curr_beat_detectors, rel_metric, rel_strategy_sets, fig_name, uParams);


%% IBIs with and without qual box plot
datasets = fieldnames(res);
if sum(strcmp(datasets, 'ppg_dalia_working'))
    primary_dataset = 'ppg_dalia_working';
else
    primary_dataset = datasets{1};
end
rel_metric = 'acc_ibi';
rel_strategy_sets = {'noQual', 'comb_noAccel'};
fig_name = 'IBI_acc_with_without_qual_boxplot';
make_boxplot_figure_two_approaches(res, primary_dataset, datasets, rel_metric, rel_strategy_sets, fig_name, uParams);

%% IBIs with comb qual box plot
datasets = fieldnames(res);
if sum(strcmp(datasets, 'ppg_dalia_working'))
    primary_dataset = 'ppg_dalia_working';
else
    primary_dataset = datasets{1};
end
rel_metric = 'acc_ibi';
rel_strategy_set = 'comb_noAccel';
fig_name = 'IBI_acc_comb_qual_boxplot';
make_boxplot_figure(res, primary_dataset, rel_metric, rel_strategy_set, fig_name, uParams);

end

function make_boxplot_figure(res, primary_dataset, curr_datasets, curr_beat_detectors, rel_metric, rel_strategy_set, fig_name, uParams)

fprintf('\n - Making boxplot figure for: %s', rel_metric);


eval(['beat_detectors = res.' primary_dataset '.noQual.num.strategy;']);
beat_detectors = curr_beat_detectors;
beat_detector_names = correct_beat_detector_names(beat_detectors);

datasets = curr_datasets;
no_y = ceil(length(datasets)/4);
figure('Position', [20,20,1000,no_y*300])
[ylims, ylab_txt] = get_ylabel_settings(rel_metric);
ftsize = 14;

for dataset_no = 1 : length(datasets)
    
    curr_dataset = datasets{dataset_no};
    
    % extract numerical results
    eval(['curr_strategy = res.' curr_dataset '.' rel_strategy_set '.num.strategy;']);
    stat_types = {'med', 'pc10', 'lq', 'uq', 'pc90'};
    for stat_type_no = 1 : length(stat_types)
        curr_stat = stat_types{stat_type_no};
        eval(['curr_metric.' curr_stat ' = res.' curr_dataset '.' rel_strategy_set '.num.' rel_metric '_' curr_stat ';']);
    end
    clear stat_type_no curr_stat stat_types
    
    % extract boxplot results for each beat detector
    for beat_detector_no = 1 : length(beat_detectors)
        curr_beat_detector = beat_detectors{beat_detector_no};
        rel_el = find(strcmp(curr_strategy, curr_beat_detector));
        data(beat_detector_no, :) = [0, curr_metric.pc10(rel_el), curr_metric.lq(rel_el), curr_metric.med(rel_el), curr_metric.uq(rel_el), curr_metric.pc90(rel_el), 100];
    end
    clear beat_detector_no rel_el curr_beat_detector curr_metric curr_strategy
    
    % create boxplot
    % inspired by: https://www.mathworks.com/matlabcentral/answers/597127-draw-a-boxplot-from-percentiles
    data = data(:, [2,2,3,3,4,4,4,5,5,6,6]);
    no_x = ceil(length(datasets)/no_y);
    if rem(length(datasets),2)==1 && dataset_no >= no_x
        temp_dataset_no = dataset_no+1;
    else
        temp_dataset_no = dataset_no;
    end
    curr_x = rem(temp_dataset_no-1, no_x)+1;
    curr_y = no_y-floor((temp_dataset_no-1)/no_x);
    subplot('Position', [0.05+(0.95*(curr_x-1)/no_x), 0.08+((curr_y-1)/no_y), 0.85/no_x, 0.6/no_y])
    boxplot(data', beat_detector_names, 'Widths',0.8, 'MedianStyle', 'target', 'Symbol', 'r')
    set(gca, 'XTickLabelRotation', 90, 'FontSize', ftsize-2)
    ylim(ylims)
    if curr_x == 1
        ylabel(ylab_txt, 'FontSize', ftsize)
    else
        set(gca, 'YTickLabel', [])
    end
    set(gca,'XGrid', 'on', 'YGrid', 'on')
    box off
    title_txt = create_title_text(curr_dataset);
    title(title_txt, 'FontSize', ftsize)
    
    % adjust whiskers to be 10th and 90th percentiles, inspired by: https://www.mathworks.com/matlabcentral/fileexchange/22526-box-plot-with-whiskers-plotted-at-fixed-percentiles
    lower_whiskers = findobj(gca,'Tag','Lower Whisker');
    lower_whisker_ends = findobj(gca,'Tag','Lower Adjacent Value');
    upper_whiskers = findobj(gca,'Tag','Upper Whisker');
    upper_whisker_ends = findobj(gca,'Tag','Upper Adjacent Value');
    for whisker_no = 1 : length(lower_whiskers)
        rel_box_no = lower_whiskers(whisker_no).XData(1);
        set(lower_whiskers(whisker_no),'YData',[data(rel_box_no,1), lower_whiskers(whisker_no).YData(2)], 'LineStyle', '-');
        set(lower_whisker_ends(whisker_no),'YData',ones(1,2)*data(rel_box_no,1));
        set(upper_whiskers(whisker_no),'YData',[upper_whiskers(whisker_no).YData(1), data(rel_box_no,11)], 'LineStyle', '-');
        set(upper_whisker_ends(whisker_no),'YData',ones(1,2)*data(rel_box_no,11));
    end
    
    % add colour to boxes, based on: https://www.mathworks.com/matlabcentral/answers/392679-how-to-fill-boxes-in-boxplot-with-different-colors
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),[0,0,1],'FaceAlpha',.5);  % greyscale: 0.5*ones(1,3)
    end
    
    % change color of median target point
    h1 = findobj(gca,'Tag','MedianInner');
    h2 = findobj(gca,'Tag','MedianInner');
    for j=1:length(h1)
        hold on
        plot(h1(j).XData, h1(j).YData, 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k')
        plot(h2(j).XData, h2(j).YData, 'ok', 'MarkerSize', 10)
        delete(h1(j))
        delete(h2(j))
    end
    
    clear data title_txt
    
end

% annotations
dim = [.165 .97 .7 .1];
str = 'Hospital Monitoring';
annotation('textbox',dim,'String',str,'FitBoxToText','on','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize', ftsize+4, 'LineStyle', 'none', 'fontweight', 'bold');
annotation('line', [.05,.415],[.985, .985])
annotation('line', [.62,.97],[.985, .985])

% dim = [.615 .97 .5 .1];
% str = 'Polysomnography Tests';
% annotation('textbox',dim,'String',str,'FitBoxToText','on','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize', ftsize+4, 'LineStyle', 'none', 'fontweight', 'bold');

dim = [.25 .64 .5 .1];
str = 'Wearable data during different emotions';
annotation('textbox',dim,'String',str,'FitBoxToText','on','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize', ftsize+4, 'LineStyle', 'none', 'fontweight', 'bold');
annotation('line', [.05,.31],[.655, .655])
annotation('line', [.69,.97],[.655, .655])

dim = [.25 .31 .5 .1];
str = 'Wearable data during activities of daily living';
annotation('textbox',dim,'String',str,'FitBoxToText','on','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize', ftsize+4, 'LineStyle', 'none', 'fontweight', 'bold');
annotation('line', [.05,.29],[.325, .325])
annotation('line', [.71,.97],[.325, .325])


do_beat_detector_annotation = 1;
if do_beat_detector_annotation
    
    dim = [.135 .0 .05 .025];
    str = 'Beat detector';
    annotation('textbox',dim,'String',str,'FitBoxToText','on','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize', ftsize, 'LineStyle', 'none');
    
    dim = [.37 .0 .05 .025];
    str = 'Beat detector';
    annotation('textbox',dim,'String',str,'FitBoxToText','on','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize', ftsize, 'LineStyle', 'none');
    
    dim = [.61 .0 .05 .025];
    str = 'Beat detector';
    annotation('textbox',dim,'String',str,'FitBoxToText','on','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize', ftsize, 'LineStyle', 'none');
    
    dim = [.85 .0 .05 .025];
    str = 'Beat detector';
    annotation('textbox',dim,'String',str,'FitBoxToText','on','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize', ftsize, 'LineStyle', 'none');
    
end

% save figure
save_figure(gcf, fig_name, uParams);
end

function make_boxplot_figure_two_approaches(res, primary_dataset, curr_datasets, rel_metric, rel_strategy_sets, fig_name, uParams)

eval(['beat_detectors = res.' primary_dataset '.noQual.num.strategy;']);

figure('Position', [20,20,1000,400])
[ylims, ylab_txt] = get_ylabel_settings(rel_metric);
ftsize = 14;
datasets = curr_datasets;

for dataset_no = 1 : length(datasets)
    
    curr_dataset = datasets{dataset_no};
    
    % extract numerical results (approaches 1 and 2)
    for approach_no = 1 : 2
        eval(['curr_strategy' num2str(approach_no) ' = res.' curr_dataset '.' rel_strategy_sets{approach_no} '.num.strategy;']);
        stat_types = {'med', 'pc10', 'lq', 'uq', 'pc90'};
        for stat_type_no = 1 : length(stat_types)
            curr_stat = stat_types{stat_type_no};
            eval(['curr_metric.' curr_stat '(:,approach_no) = res.' curr_dataset '.' rel_strategy_sets{approach_no} '.num.' rel_metric '_' curr_stat ';']);
        end
        clear stat_type_no curr_stat stat_types
    end
    
    % extract boxplot results for each beat detector
    labels = {};
    for beat_detector_no = 1 : length(beat_detectors)
        curr_beat_detector = beat_detectors{beat_detector_no};
        rel_el = find(strcmp(curr_strategy1, curr_beat_detector));
        data(2*beat_detector_no-1, :) = [0, curr_metric.pc10(rel_el,1), curr_metric.lq(rel_el,1), curr_metric.med(rel_el,1), curr_metric.uq(rel_el,1), curr_metric.pc90(rel_el,1), 100];
        labels{end+1,1} = [curr_beat_detector 'n'];
        rel_el = find(strcmp(curr_strategy2, curr_beat_detector));
        data(2*beat_detector_no, :) = [0, curr_metric.pc10(rel_el,2), curr_metric.lq(rel_el,2), curr_metric.med(rel_el,2), curr_metric.uq(rel_el,2), curr_metric.pc90(rel_el,1), 100];
        labels{end+1,1} = [curr_beat_detector 'q'];
    end
    clear beat_detector_no rel_el curr_beat_detector curr_metric curr_strategy
    
    % create boxplot
    % inspired by: https://www.mathworks.com/matlabcentral/answers/597127-draw-a-boxplot-from-percentiles
    data = data(:, [2,2,3,3,4,4,4,5,5,6,6]);
    subplot('Position', [0.05+(0.95*(dataset_no-1)/length(datasets)), 0.12, 0.85/length(datasets), 0.82])
    boxplot(data', labels, 'MedianStyle', 'target', 'Symbol', 'r')
    set(gca, 'XTickLabelRotation', 90)
    ylim(ylims)
    if dataset_no == 1
        ylabel(ylab_txt, 'FontSize', ftsize)
    else
        set(gca, 'YTickLabel', [])
    end
    set(gca,'XGrid', 'on', 'YGrid', 'on')
    box off
    title_txt = create_title_text(curr_dataset);
    title(title_txt, 'FontSize', ftsize)
    
    % adjust whiskers to be 10th and 90th percentiles, inspired by: https://www.mathworks.com/matlabcentral/fileexchange/22526-box-plot-with-whiskers-plotted-at-fixed-percentiles
    lower_whiskers = findobj(gca,'Tag','Lower Whisker');
    lower_whisker_ends = findobj(gca,'Tag','Lower Adjacent Value');
    upper_whiskers = findobj(gca,'Tag','Upper Whisker');
    upper_whisker_ends = findobj(gca,'Tag','Upper Adjacent Value');
    for whisker_no = 1 : length(lower_whiskers)
        rel_box_no = lower_whiskers(whisker_no).XData(1);
        set(lower_whiskers(whisker_no),'YData',[data(rel_box_no,1), lower_whiskers(whisker_no).YData(2)], 'LineStyle', '-');
        set(lower_whisker_ends(whisker_no),'YData',ones(1,2)*data(rel_box_no,1));
        set(upper_whiskers(whisker_no),'YData',[upper_whiskers(whisker_no).YData(1), data(rel_box_no,11)], 'LineStyle', '-');
        set(upper_whisker_ends(whisker_no),'YData',ones(1,2)*data(rel_box_no,11));
    end
    
    % add colour to boxes, based on: https://www.mathworks.com/matlabcentral/answers/392679-how-to-fill-boxes-in-boxplot-with-different-colors
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        if rem(j,2) == 0
            patch(get(h(j),'XData'),get(h(j),'YData'),[0,0,1],'FaceAlpha',.5);
        else
            patch(get(h(j),'XData'),get(h(j),'YData'),[1,0,0],'FaceAlpha',.5);
        end
    end
    
    % change color of median target point
    h1 = findobj(gca,'Tag','MedianInner');
    h2 = findobj(gca,'Tag','MedianInner');
    for j=1:length(h1)
        hold on
        plot(h1(j).XData, h1(j).YData, 'ok', 'MarkerSize', 3, 'MarkerFaceColor', 'k')
        plot(h2(j).XData, h2(j).YData, 'ok', 'MarkerSize', 6)
        delete(h1(j))
        delete(h2(j))
    end
    
    clear data title_txt
    
end

% save figure
save_figure(gcf, fig_name, uParams);
end

function make_boxplot_figure_multiple_approaches(res, primary_dataset, curr_datasets, curr_beat_detectors, rel_metric, rel_strategy_sets, fig_name, uParams)

eval(['beat_detectors = res.' primary_dataset '.noQual.num.strategy;']);
beat_detectors = curr_beat_detectors;

figure('Position', [20,20,1000,300])
[ylims, ylab_txt] = get_ylabel_settings(rel_metric);
ftsize = 14;
datasets = curr_datasets;

for dataset_no = 1 : length(datasets)
    
    curr_dataset = datasets{dataset_no};
    
    % extract numerical results (approaches 1 and 2)
    for approach_no = 1 : length(rel_strategy_sets)
        eval(['curr_strategy' num2str(approach_no) ' = res.' curr_dataset '.' rel_strategy_sets{approach_no} '.num.strategy;']);
        stat_types = {'med', 'pc10', 'lq', 'uq', 'pc90'};
        for stat_type_no = 1 : length(stat_types)
            curr_stat = stat_types{stat_type_no};
            eval(['curr_metric.' curr_stat '(:,approach_no) = res.' curr_dataset '.' rel_strategy_sets{approach_no} '.num.' rel_metric '_' curr_stat ';']);
        end
        clear stat_type_no curr_stat stat_types
    end
    
    % extract boxplot results for each beat detector
    labels = {};
    for beat_detector_no = 1 : length(beat_detectors)
        curr_beat_detector = beat_detectors{beat_detector_no};
        for approach_no = 1 : length(rel_strategy_sets)
            curr_strategy = rel_strategy_sets{approach_no};
            eval(['rel_el = find(strcmp(curr_strategy' num2str(approach_no) ', curr_beat_detector));']);
            data(length(rel_strategy_sets)*beat_detector_no-length(rel_strategy_sets)+approach_no, :) = [0, curr_metric.pc10(rel_el,approach_no), curr_metric.lq(rel_el,approach_no), curr_metric.med(rel_el,approach_no), curr_metric.uq(rel_el,approach_no), curr_metric.pc90(rel_el,approach_no), 100];
            if approach_no == 1
                % labels{end+1,1} = [curr_beat_detector '_' curr_strategy];
                labels{end+1,1} = curr_beat_detector;
            else
                labels{end+1,1} = num2str(length(labels)+1);
            end
        end
    end
    clear beat_detector_no rel_el curr_beat_detector curr_metric curr_strategy
    
    % create boxplot
    % inspired by: https://www.mathworks.com/matlabcentral/answers/597127-draw-a-boxplot-from-percentiles
    data = data(:, [2,2,3,3,4,4,4,5,5,6,6]);
    if strcmp(rel_metric, 'ppv')
        subplot('Position', [0.05+(0.95*(dataset_no-1)/length(datasets)), 0.12, 0.8/length(datasets), 0.72])
    else
        subplot('Position', [0.05+(0.95*(dataset_no-1)/length(datasets)), 0.12, 0.8/length(datasets), 0.77])
    end
    boxplot(data', labels, 'MedianStyle', 'target', 'Symbol', 'r')
    set(gca, 'XTickLabelRotation', 90)
    ylim(ylims)
    if dataset_no == 1
        ylabel(ylab_txt, 'FontSize', ftsize)
    else
        set(gca, 'YTickLabel', [])
    end
    set(gca,'XGrid', 'on', 'YGrid', 'on')
    box off
    title_txt = create_title_text(curr_dataset);
    title(title_txt, 'FontSize', ftsize)
    rel_ticks = 1:length(rel_strategy_sets):length(labels);
    set(gca, 'XTick', rel_ticks, 'XTickLabel', labels(rel_ticks))
    
    % adjust whiskers to be 10th and 90th percentiles, inspired by: https://www.mathworks.com/matlabcentral/fileexchange/22526-box-plot-with-whiskers-plotted-at-fixed-percentiles
    lower_whiskers = findobj(gca,'Tag','Lower Whisker');
    lower_whisker_ends = findobj(gca,'Tag','Lower Adjacent Value');
    upper_whiskers = findobj(gca,'Tag','Upper Whisker');
    upper_whisker_ends = findobj(gca,'Tag','Upper Adjacent Value');
    for whisker_no = 1 : length(lower_whiskers)
        rel_box_no = lower_whiskers(whisker_no).XData(1);
        set(lower_whiskers(whisker_no),'YData',[data(rel_box_no,1), lower_whiskers(whisker_no).YData(2)], 'LineStyle', '-');
        set(lower_whisker_ends(whisker_no),'YData',ones(1,2)*data(rel_box_no,1));
        set(upper_whiskers(whisker_no),'YData',[upper_whiskers(whisker_no).YData(1), data(rel_box_no,11)], 'LineStyle', '-');
        set(upper_whisker_ends(whisker_no),'YData',ones(1,2)*data(rel_box_no,11));
    end
    
    % add colour to boxes, based on: https://www.mathworks.com/matlabcentral/answers/392679-how-to-fill-boxes-in-boxplot-with-different-colors
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        h_patch(j) = patch(get(h(j),'XData'),get(h(j),'YData'),ones(1,3)*(1-(0.8*rem(j-1,length(rel_strategy_sets))/length(rel_strategy_sets))),'FaceAlpha',.5);
    end
    
    % change color of median target point
    h1 = findobj(gca,'Tag','MedianInner');
    h2 = findobj(gca,'Tag','MedianInner');
    for j=1:length(h1)
        hold on
        plot(h1(j).XData, h1(j).YData, 'ok', 'MarkerSize', 3, 'MarkerFaceColor', 'k')
        plot(h2(j).XData, h2(j).YData, 'ok', 'MarkerSize', 6)
        delete(h1(j))
        delete(h2(j))
    end
    
    % add legend
    if dataset_no == 1 && strcmp(rel_metric, 'ppv')
        rel_els = length(h_patch) - [0:length(rel_strategy_sets)-1];
        rel_h = h_patch(rel_els);
        leg_txt = create_legend_text(rel_strategy_sets);
        if strcmp(rel_metric, 'ppv') || strcmp(rel_metric, 'sens') || strcmp(rel_metric, 'acc_hr')
            %legend(rel_h, leg_txt, 'Location', 'SouthWest');
            legend(rel_h, leg_txt, 'Position', [0.5,0.86,0.1,0.2], 'Orientation', 'Horizontal', 'FontSize',12);
        else
            legend(rel_h, leg_txt, 'Location', 'NorthWest');
        end
    end
    
    clear data title_txt
    
end

% save figure
save_figure(gcf, fig_name, uParams);
end

function leg_txt = create_legend_text(rel_strategy_sets)

leg_txt = rel_strategy_sets;
leg_txt(strcmp(rel_strategy_sets, 'noQual')) = {'No quality tool'};
leg_txt(strcmp(rel_strategy_sets, 'accel')) = {'Accelerometry'};
leg_txt(strcmp(rel_strategy_sets, 'comb_noAccel')) = {'Beat detectors'};
leg_txt(strcmp(rel_strategy_sets, 'comb_Accel')) = {'Beat detectors and Accelerometry'};

end

function [ylims, ylab_txt] = get_ylabel_settings(rel_metric)
if strcmp(rel_metric, 'f1_score')
    ylims = [40,100];
    ylab_txt = 'F_1 score (%)';
elseif strcmp(rel_metric, 'ppv')
    ylims = [60,100];
    ylab_txt = 'Positive predictive value (%)';
elseif strcmp(rel_metric, 'sens')
    ylims = [0,100];
    ylab_txt = 'Sensitivity (%)';
elseif strcmp(rel_metric, 'mae_ibi')
    ylims = [0,500];
    ylab_txt = 'IBI mean absolute error (ms)';
elseif strcmp(rel_metric, 'mae_hr')
    ylims = [0,30];
    ylab_txt = 'HR mean absolute error (bpm)';
elseif strcmp(rel_metric, 'loa_hr')
    ylims = [0,50];
    ylab_txt = 'Heart rate limits of agreement (%)';
elseif strcmp(rel_metric, 'acc_hr')
    ylims = [0,100];
    ylab_txt = 'Heart rate accuracy (%)';
elseif strcmp(rel_metric, 'acc_ibi')
    ylims = [0,100];
    ylab_txt = 'Inter-beat interval accuracy (%)';
end
end

function title_txt = create_title_text(curr_dataset)

switch curr_dataset
    case 'capnobase'
        title_txt = 'CapnoBase';
    case 'ppg_dalia_sitting'
        title_txt = 'PPG-DaLiA (sitting)';
    case 'wesad_baseline'
        title_txt = 'WESAD (baseline)';
    case 'wesad_stress'
        title_txt = 'WESAD (stress)';
    case 'wesad_meditation'
        title_txt = 'WESAD (meditation)';
    case 'wesad_amusement'
        title_txt = 'WESAD (amusement)';
    case 'mimic_non_af'
        title_txt = 'MIMIC III (non-AF)';
    case 'mimic_af'
        title_txt = 'MIMIC III (AF)';
    case 'ppg_dalia_working'
        title_txt = 'PPG-DaLiA (working)';
    case 'ppg_dalia_walking'
        title_txt = 'PPG-DaLiA (walking)';
    case 'ppg_dalia_lunch_break'
        title_txt = 'PPG-DaLiA (lunch break)';
    case 'ppg_dalia_cycling'
        title_txt = 'PPG-DaLiA (cycling)';
    case 'bidmc'
        title_txt = 'BIDMC';
    case 'mimic_train_all'
        title_txt = 'MIMIC PERform (training)';
    case 'mimic_test_all'
        title_txt = 'MIMIC PERform (testing)';
    otherwise
        title_txt = strrep(curr_dataset, '_', ' ');
end

end

function save_figure(fig_handle, fig_name, uParams)

fig_path = [uParams.paths.plots_folder, fig_name];
print(fig_path, '-depsc');
close all

end