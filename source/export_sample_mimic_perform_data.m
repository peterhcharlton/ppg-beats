function export_sample_mimic_perform_data
% EXPORT_SAMPLE_MIMIC_PERFORM_DATA  Export sample data.
%   EXPORT_SAMPLE_MIMIC_PERFORM_DATA exports sample data excerpts from the
%   MIMIC Perform datasets.
%
%   # Inputs
%   None, although the script uses the following MIMIC PERform dataset files (which can be downloaded from <https://doi.org/10.5281/zenodo.6807402>):
%
%   * mimic_perform_train_n_data.mat
%   * mimic_perform_train_all_data.mat
%   * mimic_perform_af_data.mat
%   * mimic_perform_non_af_data.mat
%
%   # Outputs
%   * files : MATLAB files containing the data excerpts.
%
%   # Preparation:
%   Modify the MATLAB script by specifying the filepaths of the data files and the folder in which to save excerpts (within the setup_up function).
%
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%
%   # Author
%   Peter H. Charlton, University of Cambridge, June 2022.
%
%   # License - GPL-3.0
%      Copyright (c) 2022 Peter H. Charlton
%      This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%      This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%      You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

fprintf('\n ~~~ Exporting Sample MIMIC PERform dataset excerpts ~~~')

up = setup_up;

extract_and_save_excerpts(up);

end

function up = setup_up

fprintf('\n - Setting up parameters')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS TO MODIFY %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

up.paths.mimic_perform_af_data = '/Users/petercharlton/Documents/Data/mimic_perform_af_dataset/mimic_perform_af_data.mat';
up.paths.mimic_perform_non_af_data = '/Users/petercharlton/Documents/Data/mimic_perform_af_dataset/mimic_perform_non_af_data.mat';
up.paths.mimic_perform_train_n_data = '/Users/petercharlton/Documents/Data/mimic_perform_train_test_datasets/mimic_perform_train_n_data.mat';
up.paths.mimic_perform_train_all_data = '/Users/petercharlton/Documents/Data/mimic_perform_train_test_datasets/mimic_perform_train_all_data.mat';

up.paths.save_folder = '/Users/petercharlton/Downloads/';

%%%%%%%%%%%%%%%%%%%%%%%%   OTHER PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%

up.settings.durn_short_samples = 60; % durn of short samples in seconds
up.settings.n_subj_truncated_dataset = 10; % number of subjects to be included in the truncated dataset
up.settings.short_samples.dataset = {'mimic_perform_non_af_data', 'mimic_perform_af_data', 'mimic_perform_train_n_data', 'mimic_perform_non_af_data'};
up.settings.short_samples.subj_no = [16, 1, 1, 7];
up.settings.short_samples.start_sample = [86000, 40000, 29600, 60000];
up.settings.short_samples.name = {'normal', 'AF', 'neonate', 'noisy'};
up.settings.fs = 125; % assume that the sampling frequency is always 125 Hz (which is the case for MIMIC III).

end

function extract_and_save_excerpts(up)

fprintf('\n - Extracting and saving excerpts: ')

%% Load datasets
datasets = fieldnames(up.paths);
datasets = datasets(contains(datasets, 'mimic'));
for dataset_no = 1 : length(datasets)
    dataset_name = datasets{dataset_no};
    eval([dataset_name ' = load(up.paths.' dataset_name ');'])
end

%% Short samples
for sample_no = 1 : length(up.settings.short_samples.dataset)

    sample_name = up.settings.short_samples.name{sample_no};
    fprintf('\n    - Sample: %s', sample_name)

    dataset_name = up.settings.short_samples.dataset{sample_no};
    subj_no = up.settings.short_samples.subj_no(sample_no);
    start_sample = up.settings.short_samples.start_sample(sample_no);
    no_samples = up.settings.fs*up.settings.durn_short_samples;
    
    % extract data from this dataset
    eval(['rel_data = ' dataset_name ';']);

    % narrow down to only the required subject
    rel_data.data = rel_data.data(subj_no);

    % narrow down to only the required signal samples
    rel_els = start_sample : (start_sample + no_samples);
    signals = fieldnames(rel_data.data);
    signals = signals(~strcmp(signals, 'fix') & ~strcmp(signals, 'abp'));
    for signal_no = 1 : length(signals)
        curr_signal = signals{signal_no};
        eval(['rel_data.data.' curr_signal '.v = rel_data.data.' curr_signal '.v(rel_els);']);
    end

    % save this short sample
    filepath = [up.paths.save_folder, 'MIMIC_PERform_1_min_' sample_name];
    data = rel_data.data;
    license = rel_data.license;
    save(filepath, 'data', 'license');

end

end
