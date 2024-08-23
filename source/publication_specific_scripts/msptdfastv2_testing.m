function msptdfastv2_testing
% MSPTDFASTV2_TESTING  Used for testing of MSPTDfast v.2.
%   MSPTDFASTV2_TESTING analyses the performance
%   (execution time and F1-score) of different configurations of the MSPTD
%   PPG beat detector.
%   
%   # Inputs
%   * ppg_detect_stats.mat - results files from analysis of: PPG-DaLiA
%   (all_activities) dataset, and MIMIC PERform (testing) dataset.
%   
%   # Outputs
%   * Text summarising performance.
%   
%   # Exemplary usage
%   
%       msptdfastv2_testing      % runs the analysis
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Author
%   Peter H. Charlton, University of Cambridge, August 2024.
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

    % load data
    eval(['up.paths.stats_file = up.paths.stats_file_' curr_dataset ';']);
    res{dataset_no} = load(up.paths.stats_file);
    up.paths.raw_file = strrep(up.paths.stats_file, 'stats', 'perf');
    perf{dataset_no} = load(up.paths.raw_file);

end

% generate plot for each option
fprintf('\n\n ~~~ Summarising test results ~~~')
metrics = {'f1_score', 'perc_time'};
summarise_testing_results(res, perf, up.datasets, metrics, up);
metrics = {'sens', 'ppv'};
summarise_testing_results(res, perf, up.datasets, metrics, up);


end

function up = setup_up

up.paths.stats_file_wesad_all_activities = '/Users/petercharlton/Library/CloudStorage/GoogleDrive-peterhcharlton@gmail.com/My Drive/Work/Publications/In Preparation/2024 MSPTDfastv2/figures/benchmarking/wesad/ppg_detect_stats.mat';
up.paths.stats_file_mimic_test_all = '/Users/petercharlton/Library/CloudStorage/GoogleDrive-peterhcharlton@gmail.com/My Drive/Work/Publications/In Preparation/2024 MSPTDfastv2/figures/benchmarking/mimic/ppg_detect_stats.mat';
up.paths.plot_folder = '/Users/petercharlton/Library/CloudStorage/GoogleDrive-peterhcharlton@gmail.com/My Drive/Work/Publications/In Preparation/2024 MSPTDfastv2/figures/';

up.datasets = {'mimic_test_all', 'wesad_all_activities'};
up.dataset_names = {'MIMIC PERform (Testing)', 'WESAD (all activities)'};

up.options.name = {'algs'; ...
    };
up.options.label = {'Algorithm'; ...
    };
up.options.values = {{'MSPTD', 'MSPTDfastv2'}; ...
    };
up.options.type = {'class'; ...
    };
up.options.alg_abbr = {{'MSPTD', 'MSPTDfastv2'}; ...
    };
up.options.opt_config = {'MSPTDfastv2'; ...
    };

% % only look at final algorithm
% fields = fieldnames(up.options);
% for field_no = 1 : length(fields)
%     eval(['up.options.' fields{field_no} ' = up.options.' fields{field_no} '(end);'])
% end

end

function summarise_testing_results(res, perf, datasets, metrics, up)

% identify algorithms
algs = res{1}.res.noQual.num.strategy;

% identify reference alg
ref_alg = up.options.opt_config{1};
        
% do stats tests
for dataset_no = 1 : length(datasets)

    for metric_no = 1 : length(metrics)
        curr_metric = metrics{metric_no};

        for alg_no = 1 : length(algs)
            comparator_alg = algs{alg_no};

            % extract data
            eval(['ref_data = perf{dataset_no}.ppg_strategy_perf.raw.' ref_alg '__none.' curr_metric ';']);
            eval(['comparator_data = perf{dataset_no}.ppg_strategy_perf.raw.' comparator_alg '__none.' curr_metric ';']);

            % perform test
            p(dataset_no, alg_no, metric_no) = signrank(ref_data,comparator_data);
            med(dataset_no, alg_no, metric_no) = median(comparator_data);
            lq(dataset_no, alg_no, metric_no) = quantile(comparator_data, 0.25);
            uq(dataset_no, alg_no, metric_no) = quantile(comparator_data, 0.75);

        end

    end

end


% create table of results
alpha = 0.05;

fprintf('\n Algorithm ')
for dataset_no = 1 : length(datasets)
    fprintf('& \\multicolumn{2}{c}{%s}', up.dataset_names{dataset_no});
end
fprintf('\\\\ \\hline \n ')
for dataset_no = 1 : length(datasets)
    for metric_no = 1 : length(metrics)
        metric_label = find_metric_label(metrics{metric_no});
        fprintf('& %s ', metric_label)
    end
end
fprintf('\\\\ \\hline')

for alg_no = 1 : length(algs)
    curr_alg = algs{alg_no};
    
    curr_line = get_line_of_results(curr_alg, algs, metrics, p, med, lq, uq, alpha, datasets);
    fprintf(curr_line)
end

fprintf('\\hline \n')


end

function curr_line = get_line_of_results(curr_alg, algs, metrics, p, med, lq, uq, alpha, datasets)

% start current line for this algorithm
curr_line = sprintf('\n %s ', curr_alg);

% insert results for each dataset and each metric

for dataset_no = 1 : length(datasets)
    for metric_no = 1 : length(metrics)

        curr_metric = metrics{metric_no};
        alg_no = find(strcmp(algs, curr_alg));
        rel_p = p(dataset_no, alg_no,metric_no);
        if rel_p < alpha
            rel_p_txt = '*';
        else
            rel_p_txt = '';
        end
        if strcmp(curr_metric, 'f1_score') || strcmp(curr_metric, 'sens') || strcmp(curr_metric, 'ppv')
            no_dp = 1;
        else
            no_dp = 4;
        end
        formatting_txt = sprintf('%%.%df (%%.%df-%%.%df)', no_dp, no_dp, no_dp);
        rel_val = sprintf(formatting_txt, med(dataset_no,alg_no,metric_no), lq(dataset_no,alg_no,metric_no), uq(dataset_no,alg_no,metric_no));

        % add to current line
        curr_line = sprintf('%s & %s %s ', curr_line, rel_val, rel_p_txt);
    end

end

% end current line
curr_line = [curr_line, '\\\\'];

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

function metric_label = find_metric_label(curr_metric)

switch curr_metric
    case 'f1_score'
        metric_label = '\fscore~(\%)';
    case 'perc_time'
        metric_label = 'Execution time (\%)';
    case 'sens'
        metric_label = 'Sensitivity (\%)';
    case 'ppv'
        metric_label = 'PPV (\%)';


end

end