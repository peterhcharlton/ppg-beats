function collate_mimic_perform_train_test_datasets
% COLLATE_MIMIC_PERFORM_TRAIN_TEST_DATASETS  Collates training and testing datasets.
%   COLLATE_MIMIC_PERFORM_TRAIN_TEST_DATASETS downloads and collates data from the MIMIC-III
%   Waveform Database, for PPG beat detector performance evaluation.
%
%   # Inputs
%   * none (although see 'setup_up' for the parameters to be set)
%
%   # Working Files
%   The script downloads several files from PhysioNet, which are saved locally.
%
%   # Outputs
%   * files : MATLAB files containing the required data for PPG beat detector evaluation: 'mimic_train_a_data.mat', 'mimic_train_n_data.mat', 'mimic_train_all_data.mat', 'mimic_test_a_data.mat', 'mimic_test_n_data.mat', and 'mimic_test_all_data.mat'.
%
%   # Preparation
%   Modify the MATLAB script by inserting the 'up.paths.root_folder' and 'up.paths.save_folder' into the setup_up function.
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

fprintf('\n ~~~ Downloading and collating MIMIC PERform Train and Test Datasets ~~~\n')

up = setup_up;

download_data(up);

end

function up = setup_up

fprintf('\n - Setting up parameters:')

% The following root folder contains subfolders for each subject (name
% 'S#' where # is the subject number)
up.paths.local.root_folder = '/Users/petercharlton/Documents/Data/mimic_perform_train_test_datasets/';
if ~exist(up.paths.local.root_folder, 'dir'), warning('Specified folder does not exist.\n Please adjust ''up.paths.local.root_folder'' in ''setup_up'''); end
fprintf('\n    - Folder in which to save data:\n    - %s', up.paths.local.root_folder)

% Location of files online
up.paths.web.database_folder = 'mimic3wdb/1.0/';
up.paths.web.root_folder = ['https://physionet.org/files/', up.paths.web.database_folder];
fprintf('\n    - Website from which to download data:\n    - %s', up.paths.web.root_folder)
up.paths.filenames.records_numerics_file = 'RECORDS-numerics';
up.paths.web.records_numerics_file = [up.paths.web.root_folder, up.paths.filenames.records_numerics_file];
up.paths.filenames.records_waveforms_file = 'RECORDS-waveforms';
up.paths.web.records_waveforms_file = [up.paths.web.root_folder, up.paths.filenames.records_waveforms_file];
up.paths.filenames.records_adults_file = 'RECORDS-adults';
up.paths.web.records_adults_file = [up.paths.web.root_folder, up.paths.filenames.records_adults_file];
up.paths.filenames.records_neonates_file = 'RECORDS-neonates';
up.paths.web.records_neonates_file = [up.paths.web.root_folder, up.paths.filenames.records_neonates_file];

% Filepaths at which to save
up.paths.local.temp_folder = [up.paths.local.root_folder, 'downloaded_files', filesep];
if ~exist(up.paths.local.temp_folder, 'dir')
    fprintf('\n - Making directory in which to store downloaded data at:');
    fprintf('\n    - %s', up.paths.local.temp_folder)
    mkdir(up.paths.local.temp_folder);
end
up.paths.local.records_waveforms_file = [up.paths.local.temp_folder, up.paths.filenames.records_waveforms_file];
up.paths.local.records_numerics_file = [up.paths.local.temp_folder, up.paths.filenames.records_numerics_file];
up.paths.local.records_neonates_file = [up.paths.local.temp_folder, up.paths.filenames.records_neonates_file];
up.paths.local.records_adults_file = [up.paths.local.temp_folder, up.paths.filenames.records_adults_file];

% Settings
up.settings.req_durn = 1*60*10; % minimum duration of recording (in secs)
up.settings.no_subjs_per_group_per_dataset = 2; %100; %%%%% CHANGE
up.settings.ecg_labels = {'II', 'I', 'V', 'III', 'MCL1'};
up.settings.datasets = {'train', 'test'};
up.settings.no_subjs_per_group = up.settings.no_subjs_per_group_per_dataset * length(up.settings.datasets);
fprintf('\n - Collating training and testing datasets, each of \n   - %.1d minute recordings\n   - %d subjects (%d adults and %d neonates)', up.settings.req_durn/60, up.settings.no_subjs_per_group_per_dataset*2, up.settings.no_subjs_per_group_per_dataset, up.settings.no_subjs_per_group_per_dataset);

% - check WFDB toolbox installation
check_WFDB_installation;

end

function check_WFDB_installation

% Check whether the 'wrsamp' function is available
if ~exist('wrsamp.m', 'file')
    error('This code requires the Waveform Database Software Package (WFDB) for MATLAB and Octave, available from: https://doi.org/10.13026/6zcz-e163 . Please install it and check it is visible in the Matlab path.');
else
    fprintf('\n    - Detected that the Waveform Database Software Package (WFDB) for MATLAB and Octave is installed.')
end

end

function download_data(up)

%% Identify records which meet the requirements
% To be included, records must contain:
% - at least 10 mins continuously
% - of the required variables: (i) Pleth, (ii) ECG, and (iii) a respiration signal.

stay_list_filepath = [up.paths.local.temp_folder, 'stay_list.mat'];
if ~exist(stay_list_filepath, 'file')
    stay_list = define_stays_for_download(up);
    save(stay_list_filepath, 'stay_list');
else
    load(stay_list_filepath)
end

%% Download data for each subject which meets the requirements

fprintf('\n - Downloading data for each subject which meets the requirements')
% cycle through each subject
inc_subj_no = 0;
groups = {};
for subj_no = 1 : length(stay_list)

    % check to see whether there are already enough of this subject's group
    curr_group = stay_list(subj_no).group;
    no_of_this_group = sum(strcmp(groups, curr_group));
    if no_of_this_group >= up.settings.no_subjs_per_group
        continue
    end

    curr_subj = stay_list(subj_no).stayID;
    fprintf('\n   - Extracting data for subject no. %d (%s)', subj_no, curr_subj)

    % check whether this subject should be skipped (subjects which were manually identified as having unsuitable data)
    if sum(strcmp(curr_subj, {'3004104','3091963','3110548','3024047','3105058', '3001551', '3002751', '3003405', '3005850', '3006240', '3013467', '3015738', '3021205', '3024873', '3025163', '3026618', '3031584', '3034102', '3039719', '3061958', '3085535', '3098373', '3124490', '3130412', '3136247', '3142279', '3027856', '3007696'}))
        fprintf('\n      - unsuitable data')
        continue
    end

    % identify recording files for this subject
    curr_rec = stay_list(subj_no).file;
    rec_file_info = get_rec_file_info(curr_rec, up);
    if ~isstruct(rec_file_info)
        fprintf('\n      - couldn''t identify recording files in recording %s', curr_rec)
        continue
    end

    % identify recording files within the specified period
    rec_file_info = get_rec_files_in_period(rec_file_info);

    % MIMIC Waveform database subfolder for this subject
    curr_rec_subfolder = curr_rec(1:2);

    % see whether each file for this recording contains the required signals
    for file_no = 1 : min([25, length(rec_file_info.rec_deb)])

        fprintf('; file %d', file_no);
        inc_curr_file_data = false;
        curr_file = rec_file_info.rec_name{file_no};

        % - check whether this file is of a suitable duration
        file_durn = (rec_file_info.rec_no_samps(file_no)-1)/rec_file_info.fs; % in secs
        if file_durn < up.settings.req_durn
            fprintf(' (duration too short)');
            continue
        end
        if file_durn > 3*up.settings.req_durn
            fprintf(' (duration too long)');
            continue
        end

        % - find out what signals this file contains
        vars = extract_details_of_signals(curr_subj, curr_rec, curr_rec_subfolder, curr_file, up);

        % - check whether this file has the desired signals
        if ~sum(contains(vars, {'II', 'I', 'V', 'III', 'MCL1'})) || ~sum(contains(vars, 'PLETH')) || ~sum(contains(vars, 'RESP'))
            fprintf(' (didn''t have the signals)');
            continue
        end

        % - download data for this file (to downloaded files folder)
        fprintf(' - downloading,')
        cd(up.paths.local.temp_folder)
        filename = [curr_file, 'm'];
        if ~exist([filename, '.mat'], 'file')
            wfdb2mat([up.paths.web.database_folder, curr_rec_subfolder, '/', curr_rec(1:7), '/', curr_file]);
        end

        % - extract information on this file
        try
            [tm,signal,Fs,labels]=rdmat(filename);
        catch
            fprintf(' (couldn''t get info on file)')
            continue
        end
        sig_names = extractfield(labels, 'Description');

        % - extract required data from this file
        curr_file_data = extract_data_from_file(sig_names, signal, Fs, curr_subj, curr_rec, curr_group, filename, up);
        fprintf(' %d mins', round(length(curr_file_data.ppg.v)/(60*curr_file_data.ppg.fs)) );

        % - check whether all sigs were of the same duration
        waves = {'ppg', 'ekg', 'imp'};
        durns = nan(size(waves));
        for wave_no = 1 : length(waves)
            eval(['durns(wave_no) = length(curr_file_data.' waves{wave_no} '.v);']);
        end
        if length(unique(durns)) > 1
            fprintf(' (signals not same duration)')
            continue
        end

        % - Refine duration of data to give required duration
        curr_file_data = refine_duration(curr_file_data, up);
        if isstr(curr_file_data)
            fprintf(' %s', curr_file_data);
            continue
        else
            fprintf(' reduced to %d mins', round(length(curr_file_data.ppg.v)/(60*curr_file_data.ppg.fs)) );
        end

        % - check to see if signals are suitable
        inc_curr_file_data = check_signals(curr_file_data, up);
        if ~inc_curr_file_data
            fprintf(' (flat line)');
            continue
        end

        % - if it gets this far, then the signals are suitable, so don't search any more files for this subject
        break

    end

    % - move to next subj if there were no suitable data
    if ~inc_curr_file_data
        continue
    end

    % store the current file data
    inc_subj_no = inc_subj_no + 1;
    data(inc_subj_no) = curr_file_data;

    % calculate number of subjects of each group
    groups = cell(length(data),1);
    for s = 1 : length(groups)
        groups{s} = data(s).fix.group;
    end
    fprintf('\n     - Current number of subjects in the dataset: %d', length(data))
    unique_groups = unique(groups);
    completed_extract = true;
    % check how many subjects are in each group
    for s = 1 : length(unique_groups)
        curr_group = unique_groups{s};
        no_in_curr_group = sum(strcmp(groups, curr_group));
        fprintf('\n        - %s: %d', curr_group, no_in_curr_group);
        % check to see whether there are enough subjects of this group in the dataset
        if no_in_curr_group ~= up.settings.no_subjs_per_group
            completed_extract = false;
        end
    end
    % - skip remaining subjects if there are enough in the dataset
    if completed_extract
        break
    end

end
clear groups

%% Save data

% add source details
source.date_of_conversion = datetime('now');
source.matlab_conversion_script = [mfilename('fullpath'), '.m'];
source.raw_data_path = up.paths.local.temp_folder;

% add license details
license.readme = 'Use the following command to view the license: ''fprintf(license.details)''';
license.details = ['This dataset is licensed under the Open Data Commons Open Database License v1.0 (ODbL 1.0 license).' ...
    '\nFurther details are available at: https://opendatacommons.org/licenses/odbl/summary/' ...
    '\n\nThis dataset is derived from the MIMIC III Waveform Database:' ...
    '\nMoody, B., Moody, G., Villarroel, M., Clifford, G. D., & Silva, I. (2020). MIMIC-III Waveform Database (version 1.0). PhysioNet. https://doi.org/10.13026/c2607m' ...
    '\n\nThe MIMIC III Waveform Database is licensed under the ODbL 1.0 license.' ...
    '\n\nThe MIMIC-III database is described in the following publication:' ...
    '\nJohnson, A. E. W., Pollard, T. J., Shen, L., Lehman, L. H., Feng, M., Ghassemi, M., Moody, B., Szolovits, P., Celi, L. A., & Mark, R. G. (2016). MIMIC-III, a freely accessible critical care database. Scientific Data, 3, 160035. https://doi.org/10.1038/sdata.2016.35' ...
    '\n\nIt is available on PhysioNet, https://physionet.org/ :' ...
    '\nGoldberger, A., Amaral, L., Glass, L., Hausdorff, J., Ivanov, P. C., Mark, R., ... & Stanley, H. E. (2000). PhysioBank, PhysioToolkit, and PhysioNet: Components of a new research resource for complex physiologic signals. Circulation [Online]. 101 (23), pp. e215–e220. https://doi.org/10.1161/01.cir.101.23.e215'];

% identify subjects of different groups
groups = cell(length(data),1);
for s = 1 : length(groups)
    groups{s} = data(s).fix.group;
end

% save data for subjects of different groups in turn
fprintf('\n   - Saving data: ')
orig_data = data; clear data
% - cycle through groups and datasets
unique_groups = unique(groups);
unique_groups(end+1) = {'all'}; % an extra group combining both adult and neonate data
for s = 1 : length(unique_groups)
    curr_group = unique_groups{s};
    
    % - identify data in this group
    if ~strcmp(curr_group, 'all')
        rel_subjs = strcmp(groups, curr_group);
    else
        rel_subjs = true(length(groups),1);
        rel_subjs_a = find(strcmp(groups, 'a'));
        rel_subjs_n = find(strcmp(groups, 'n'));
    end
    rel_subjs = find(rel_subjs);

    % - save as separate datasets
    no_in_this_group_and_dataset = floor(length(rel_subjs)/length(up.settings.datasets));
    for dataset_no = 1 : length(up.settings.datasets)
        curr_dataset = up.settings.datasets{dataset_no};
        if ~strcmp(curr_group, 'all')
            curr_rel_subjs = rel_subjs((dataset_no-1)*no_in_this_group_and_dataset+1:dataset_no*no_in_this_group_and_dataset);
        else
            no_in_this_group_and_dataset_and_cohort = floor(no_in_this_group_and_dataset/2);
            curr_rel_subjs = [rel_subjs_a((dataset_no-1)*no_in_this_group_and_dataset_and_cohort+1:dataset_no*no_in_this_group_and_dataset_and_cohort); ...
                rel_subjs_n((dataset_no-1)*no_in_this_group_and_dataset_and_cohort+1:dataset_no*no_in_this_group_and_dataset_and_cohort)];
        end
        data = orig_data(curr_rel_subjs);
        for new_data_s = 1:length(data), new_groups{new_data_s,1} = data(new_data_s).fix.group; end
        dataset_name = ['mimic_' curr_dataset, '_' curr_group];
        fprintf('\n       - Dataset %s: %d adults; %d neonates', dataset_name, sum(strcmp(new_groups, 'a')), sum(strcmp(new_groups, 'n')))
        save([up.paths.local.root_folder, dataset_name '_data'], 'data', 'source', 'license')
        clear new_groups curr_dataset curr_rel_subjs no_in_this_group_and_dataset_and_cohort new_data_s dataset_name
    end
    clear rel_subjs curr_group no_in_this_group_and_dataset
end

%% Closing messages

fprintf('\n\n - NB: the raw data files also contains additional variables which have not been imported');

fprintf('\n\n ~~~ DONE ~~~')

end

function curr_file_data = extract_data_from_file(sig_names, signal, Fs, curr_subj, curr_rec, curr_group, filename, up)

%    - PPG
rel_col = find_rel_col(sig_names, 'PLETH');
curr_file_data.ppg.v = signal(:,rel_col);
curr_file_data.ppg.fs = Fs;
curr_file_data.ppg.method = 'PPG recorded using bedside monitor'; % From description at: https://physionet.org/content/mimicdb/1.0.0/

%    - ECG
for ecg_label_no = 1 : length(up.settings.ecg_labels)
    curr_ecg_label = up.settings.ecg_labels{ecg_label_no};
    rel_col = find_rel_col(sig_names, curr_ecg_label);
    if sum(rel_col)
        curr_file_data.ekg.v = signal(:,rel_col);
        curr_file_data.ekg.fs = Fs;
        curr_file_data.ekg.method = ['ECG recorded using bedside monitor, lead ', curr_ecg_label];
        break
    end
end

%    - RESP
rel_col = find_rel_col(sig_names, 'RESP');
if ~isempty(rel_col)
    curr_file_data.imp.v = signal(:,rel_col);
    curr_file_data.imp.fs = Fs;
    curr_file_data.imp.method = 'Impedance pneumography respiratory signal recorded at the chest using bedside monitor';
end

% - fixed data
curr_file_data.fix.subj_id = curr_subj;
curr_file_data.fix.rec_id = curr_rec;
curr_file_data.fix.file = filename;
curr_file_data.fix.group = curr_group;

end

function vars = extract_details_of_signals(curr_subj, curr_rec, curr_rec_subfolder, curr_file, up)
use_wfdb_toolbox = 0;
if use_wfdb_toolbox
    [siginfo,~,~] = wfdbdesc(['mimic3wdb/matched/', curr_rec_subfolder, '/', curr_rec(1:7), '/', curr_file]);
    vars = extractfield(siginfo, 'Description');
else
    rec_header_path = [up.paths.web.root_folder, curr_rec_subfolder, '/', curr_subj, '/', curr_file '.hea'];
    vars = get_vars_in_rec(rec_header_path);
end

end

function curr_file_data = refine_duration(curr_file_data, up)

sigs = fieldnames(curr_file_data);
sigs = sigs(~strcmp(sigs, 'fix'));

% - identify els which are a nan for at least one of the three signals
nan_els = false(size(curr_file_data.ekg.v));
for sig_no = 1 : length(sigs)
    eval(['rel_data = curr_file_data.' sigs{sig_no} ';']);
    sig_len = length(rel_data.v);
    if ~isempty(rel_data)
        nan_els(isnan(rel_data.v)) = true;
    end
end

% - identify a segment of non-nans:
nan_els = find(nan_els);
segs.deb = [1; nan_els+1];
segs.fin = [nan_els; sig_len];
segs.len = segs.fin-segs.deb;
req_samps = up.settings.req_durn*125;
rel_seg = find(segs.len> req_samps,1);
if isempty(rel_seg)
    curr_file_data = ' (insufficient duration without nans)';
    return
end

% - truncate data to be of required duration
for sig_no = 1 : length(sigs)
    eval(['rel_data = curr_file_data.' sigs{sig_no} ';']);
    if ~isempty(rel_data)
        eval(['curr_file_data.' sigs{sig_no} '.v = curr_file_data.' sigs{sig_no} '.v(segs.deb(rel_seg):segs.deb(rel_seg)+req_samps);']);
    end
end

end

function inc_curr_file_data = check_signals(curr_file_data, up)

sigs = fieldnames(curr_file_data);
sigs = sigs(~strcmp(sigs, 'fix'));

% - check for flat line
inc_curr_file_data = true;
for sig_no = 1 : length(sigs)
    eval(['curr_data = curr_file_data.' sigs{sig_no} '.v;']);
    eval(['curr_fs = curr_file_data.' sigs{sig_no} '.fs;']);
    
    % Adapted from: https://uk.mathworks.com/matlabcentral/answers/382011-how-to-count-the-number-of-consecutive-identical-elements-in-both-the-directions-in-a-binary-vecto#answer_304569
    Q = find([false; diff(curr_data)~=0]);
    I = zeros(length(curr_data),1) ;
    I(Q) = 1;
    I = cumsum(I);
    N = diff([1; Q; length(curr_data)+1]);
    no_consecutive = N(I+1);
    
    if strcmp(sigs{sig_no}, 'imp')
        thresh_durn = 10;
    else
        thresh_durn = 5;
    end

    if sum(curr_data == mode(curr_data)) > (0.2*length(curr_data))
        inc_curr_file_data = false;
        continue
    end

    if max(no_consecutive) > (thresh_durn*curr_fs)
        inc_curr_file_data = false;
    end
end

end

function rel_col = find_rel_col(sig_names, rel_sig)

rel_col = strcmp(sig_names, rel_sig);

if sum(rel_col)>1
    rel_col(find(rel_col,1)+1:end) = 0;
end

end

function vars = get_vars_in_rec(rec_header_path)

sig_data = webread(rec_header_path);
sig_data = splitlines(sig_data);
vars = {};

% First line contains the time and durn of recording
while strcmp(sig_data{1}(end), ' '), sig_data{1} = sig_data{1}(1:end-1); end

% Extract variable names from remaining lines
for line_no = 2 : length(sig_data)-1
    temp = strfind(sig_data{line_no}, ' ');
    if strcmp(sig_data{line_no}(end), ' ')
        vars{end+1,1} = sig_data{line_no}(temp(end-1)+1 : temp(end)-1);
    else
        vars{end+1,1} = sig_data{line_no}(temp(end)+1 : end);
    end
end


end

function rec_file_info = get_rec_file_info(curr_rec, up)

%% import header file for stay
curr_rec_subfolder = curr_rec(1:2);
curr_rec_hea_path = [up.paths.web.root_folder, curr_rec_subfolder, '/', curr_rec(1:7), '/', curr_rec, '.hea'];
try
    rec_data = webread(curr_rec_hea_path);
catch
    rec_file_info = 'couldn''t download file';
    return
end
rec_data = splitlines(rec_data);
temp = strfind(rec_data{1}, ' ');
while strcmp(rec_data{1}(end), ' '), rec_data{1} = rec_data{1}(1:end-1); end

%% Extract Information

% header line
line_no = 1;
temp = strsplit(rec_data{line_no}, ' ');
temp2 = strfind(temp{line_no}, '/');
if ~isempty(temp2)
    rec_file_info.stay_name = temp{line_no}(1:temp2-1);
    rec_file_info.no_segs = temp{line_no}(temp2+1:end);
else
    rec_file_info.stay_name = temp{line_no};
end
rec_file_info.no_sigs = str2double(temp{2});
rec_file_info.fs = str2double(temp{3});
rec_file_info.no_samps = str2double(temp{4});
rec_file_info.deb = datenum(temp{5}, 'HH:MM:SS.FFF');
rec_file_info.fin = rec_file_info.deb + (rec_file_info.no_samps/(60*60*24*rec_file_info.fs));
clear temp temp2

% remaining lines
no_samps_passed = 0;
rec_file_info.rec_name = {};
[rec_file_info.rec_deb, rec_file_info.rec_fin] = deal([]);
rec_no = 0; prev_seg_blank = 0;
while line_no < length(rec_data)
    line_no = line_no+1;

    curr_line = rec_data{line_no};

    % skip an empty line
    if isempty(curr_line)
        continue
    end

    % skip a comment line
    if strcmp(curr_line(1), '#')
        continue
    end

    % skip a layout line
    if contains(curr_line, 'layout')
        continue
    end

    % blank segments (i.e. no data in this segment)
    if strcmp(curr_line(1), '~')
        no_samps_passed = no_samps_passed + str2double(curr_line(3:end));
        prev_seg_blank = 1;
        if rec_no > 0
            rec_file_info.next_blank(rec_no,1) = 1;
        end
        continue
    end

    % recordings (i.e. segments containing data)
    rec_no = rec_no+1;
    temp = strfind(curr_line, ' ');
    rec_file_info.rec_name{rec_no,1} = curr_line(1:temp(1)-1);
    rec_file_info.rec_deb(rec_no,1) = rec_file_info.deb + datenum(no_samps_passed/(60*60*24*rec_file_info.fs));
    rec_file_info.rec_no_samps(rec_no,1) = str2double(curr_line(temp(1)+1:end));
    no_samps_passed = no_samps_passed + rec_file_info.rec_no_samps(rec_no,1);
    rec_file_info.rec_fin(rec_no,1) = rec_file_info.deb + datenum(no_samps_passed/(60*60*24*rec_file_info.fs));
    rec_file_info.prev_blank(rec_no,1) = prev_seg_blank;
    rec_file_info.next_blank(rec_no,1) = 0;

    prev_seg_blank = 0;
end

end

function rec_file_info = get_rec_files_in_period(rec_file_info)

%% Calculate start and end of period
% - start time of relevant period for this subject
period_deb = rec_file_info.deb;
% - end time of relevant period for this subject
period_fin = inf;

%% Identify recordings which are completely within this period
rec_file_info.in_period = false(size(rec_file_info.rec_deb));
for s = 1 : length(rec_file_info.rec_deb)

    if rec_file_info.rec_deb(s) >= period_deb && rec_file_info.rec_fin(s) <= period_fin
        rec_file_info.in_period(s) = 1;
    end

end

end

function records = define_stays_for_download(up)

fprintf('\n - Defining stays for download')

% This function defines which MIMIC III records should be considered for inclusion.

% Download required records files
file_types = {'neonates', 'adults', 'numerics', 'waveforms'};
for file_no = 1 : length(file_types)
    curr_file_type = file_types{file_no};
    eval(['local_filepath = up.paths.local.records_' curr_file_type '_file;']);
    if ~exist(local_filepath, 'file')
        eval(['web_filepath = up.paths.web.records_' curr_file_type '_file;']);
        websave(local_filepath, web_filepath);
    end
end

% Extract subject IDs for each group from their respective records files
fprintf('\n   - Extracting relevant subject IDs from records files');
records = struct;
counter = 0;
for group = {'neonates', 'adults'}
    curr_group = group{1,1};
    eval(['curr_filepath = up.paths.local.records_' curr_group '_file;']);
    fid = fopen(curr_filepath);
    while ~feof(fid)
        counter = counter + 1;
        records(counter).path = fgetl(fid);
        records(counter).subj_id = records(counter).path(5:end-1);
        records(counter).group = curr_group(1);
    end
    fclose(fid);
    clear fid curr_filepath curr_group
end

% Narrow down to only subject IDs who have waveforms and numerics
fprintf('\n   - Narrowing down to only include records with waveforms and numerics')
for group = {'waveforms', 'numerics'}
    curr_group = group{1,1};

    % - extract list of subjects with this data type
    eval(['curr_filepath = up.paths.local.records_' curr_group '_file;']);
    curr_records = cell(0);
    curr_subj_ids = cell(0);
    fid = fopen(curr_filepath);
    while ~feof(fid)
        curr_records{end+1,1} = fgetl(fid);
        curr_subj_ids{end+1,1} = curr_records{end}(5:10);
    end
    fclose(fid);
    clear fid curr_filepath curr_line

    if strcmp(curr_group, 'numerics')
        numerics_records = curr_records;
    end
    clear curr_records

    % - narrow down list of records to only those with this data type
    [~,records_to_keep,~] = intersect(extractfield(records, 'subj_id'), curr_subj_ids);
    records = records(records_to_keep);
    clear records_to_keep curr_subj_ids

end

% shorten the list of records to only the first 10,000, so that the next step doesn't take too long.
records = records(1:10000); % used to be 2,500

% Identify corresponding records (just taking the first one for each subject)
fprintf('\n   - Identifying the first record for each subject')
% - use the list of numerics records produced in the previous step: numerics_records
% - identify first numerics record for each subj
for record_no = 1 : length(records)
    curr_subj_id = records(record_no).subj_id;
    rel_numerics_record_el = find(contains(numerics_records, curr_subj_id),1);
    if ~isempty(rel_numerics_record_el)
        rel_numerics_file = numerics_records{rel_numerics_record_el};
        temp = strfind(rel_numerics_file, '/');
        records(record_no).file = rel_numerics_file(temp(2)+1:end-1); % excluded the 'n' from the end
    end
end
% - exclude records without a file
els_to_exclude = cellfun(@isempty, extractfield(records, 'file'));
records = records(~els_to_exclude);

% Extract record ID
for record_no = 1 : length(records)
    records(record_no).stayID = records(record_no).file(1:7);
end

end
