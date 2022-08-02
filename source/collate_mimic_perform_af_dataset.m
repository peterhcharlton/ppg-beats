function collate_mimic_perform_af_dataset
% COLLATE_MIMIC_PERFORM_AF_DATASET  Collate AF dataset.
%   COLLATE_MIMIC_PERFORM_AF_DATASET downloads and collates data from the MIMIC-III Waveform
%   Database Matched Subset during Atrial Fibrillation (AF) and sinus rhythm.
%
%   # Inputs
%   * none (although see 'setup_up' for the parameters to be set)
%
%   # Working Files
%   The script downloads several files from PhysioNet, which are saved locally.
%
%   # Outputs
%   * files : Two MATLAB files containing the required data for PPG beat detector evaluation: 'mimic_af_data.mat', and 'mimic_non_af_data.mat'.
%
%   # Preparation:
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

fprintf('\n ~~~ Downloading and collating MIMIC III matched waveform database excerpt ~~~')

up = setup_up;

download_data(up);

end

function up = setup_up

fprintf('\n - Setting up parameters')

% The following root folder contains subfolders for each subject (name
% 'S#' where # is the subject number)
up.paths.local.root_folder = '/Users/petercharlton/Documents/Data/mimiciii_ppg_af_beat_detection/';
if ~exist(up.paths.local.root_folder, 'dir'), warning('Specified folder does not exist.\n Please adjust ''up.paths.local.root_folder'' in ''setup_up'''); end

% Location of files online
up.paths.web.root_folder = 'https://physionet.org/files/mimic3wdb-matched/1.0/';

% Filepaths at which to save
up.paths.local.temp_folder = [up.paths.local.root_folder, 'downloaded_files', filesep];
if ~exist(up.paths.local.temp_folder, 'dir')
    fprintf('\n - Making directory in which to store downloaded data');
    mkdir(up.paths.local.temp_folder);
end

% Settings
up.settings.req_durn = 1*20*60; % minimum duration of recording (in secs)

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

fprintf('\n - Downloading data')

%% Identify records which meet the requirements
% To be included, records must contain:
% - at least 1 hour continuously
% - of the required variables: (i) Pleth, and (ii) an ECG signal.

stay_list = define_stays_for_download;

%% Download data for each subject which meets the requirements

% cycle through each subject
inc_subj_no = 0;
for subj_no = 1 : length(stay_list.stayID)

    % extract current subject
    curr_subj = stay_list.stayID{subj_no};

    fprintf('\n   - Extracting data for subject %s:', curr_subj)

    % identify recording files for this subject
    curr_rec = stay_list.file{subj_no};
    rec_file_info = get_rec_file_info(curr_rec, up);

    % identify recording files within the specified period
    rec_file_info = get_rec_files_in_period(rec_file_info, stay_list.onset_time{subj_no}, stay_list.offset_time{subj_no});

    curr_rec_subfolder = curr_rec(1:3);

    % see whether each file for this recording contains the required signals
    possible_files = {}; curr_durn = 0; at_least_one_file_had_signals = false;
    for file_no = 1 : length(rec_file_info.rec_deb)

        if strcmp(curr_subj, 'p092846') && file_no < 4
            continue
        end

        fprintf(', file %d', file_no);

        curr_file = rec_file_info.rec_name{file_no};

        % - extract information on this file
        use_wfdb_toolbox = 0;
        if use_wfdb_toolbox
            [siginfo,~,~] = wfdbdesc(['mimic3wdb/matched/', curr_rec_subfolder, '/', curr_rec(1:7), '/', curr_file]);
            vars = extractfield(siginfo, 'Description');
        else
            rec_header_path = [up.paths.web.root_folder, curr_rec_subfolder, '/', curr_subj, '/', curr_file '.hea'];
            vars = get_vars_in_rec(rec_header_path);
        end

        % - check whether this file has the desired signals
        if ~sum(contains(vars, {'II', 'I', 'V', 'III', 'MCL1'})) || ~sum(contains(vars, 'PLETH'))
            fprintf(' (didn''t have the signals)');
            possible_files = {};
            continue
        end
        at_least_one_file_had_signals = true;

%         if ~use_wfdb_toolbox
%             [siginfo,~,~] = wfdbdesc(['mimic3wdb/matched/', curr_rec_subfolder, '/', curr_rec(1:7), '/', curr_file]);
%             vars = extractfield(siginfo, 'Description');
%         end
        t_deb = rec_file_info.rec_deb(file_no);

        % - see whether this file has the desired start time (if not the first file to be found)
        if ~isempty(possible_files)
            % calculate the difference between this start time and the desired start time
            t_offset = (next_t_deb - t_deb)*24*60*60; % in secs
            % if the difference isn't one sample then skip this file
            if round(1/t_offset) ~= rec_file_info.fs
                fprintf(' (didn''t have the required start time)');
                possible_files = {};
                continue
            end
        else
            % - note the start of this file as a potential overall start time
            t_overall_deb = t_deb;
        end

        % - if it does have the desired signals, then note this file as a possible file
        possible_files(end+1,1) = {curr_file};

        % - see whether enough data has accumulated over the files
        curr_durn = round(60*60*24*(rec_file_info.rec_fin(file_no)-t_overall_deb));
        if curr_durn >= up.settings.req_durn
            % have got enough data from this recording, so skip the rest of the files
            break
        end

        % - calculate the desired next start time
        next_t_deb = rec_file_info.rec_fin(file_no) + ((1/rec_file_info.fs)/(24*60*60));

    end

    % Check to see whether there was enough data in this recording
    if curr_durn < up.settings.req_durn
        % - if not then skip this recording
        if at_least_one_file_had_signals
            fprintf('\n      - not enough data in recording %s (only %.1f mins)', curr_rec, curr_durn/60)
        else
            fprintf('\n      - recording %s didn''t have the required signals', curr_rec)
        end
        continue
    end

    % Download data
    fprintf('\n     - Downloading relevant data for recording %s:', curr_rec)
    inc_subj_no = inc_subj_no+1;
    % - store in downloaded files folder
    cd(up.paths.local.temp_folder)

    % - download each signal file in turn
    no_samps_required = up.settings.req_durn*rec_file_info.fs + 1;
    no_samps_downloaded = 0;
    for file_no = 1 : length(possible_files)

        if subj_no == 35 && file_no == 1 || subj_no == 22 && file_no < 5 
            continue
        end

        fprintf(' file %d,', file_no)

        curr_file = possible_files{file_no};

        % - find out how many samples to download from this file
        rel_row = strcmp(rec_file_info.rec_name, curr_file);
        no_samps_in_file = rec_file_info.rec_no_samps(rel_row);
        remaining_samps_required = no_samps_required-no_samps_downloaded;
        no_samps_to_download = min([no_samps_in_file, remaining_samps_required]);

        % - download this file
        filename = [curr_file, 'm'];
        if ~exist([filename, '.mat'], 'file')
            % I couldn't get it to work when specifying the number of samples:
            specify_no_samps = 0;
            if specify_no_samps
                wfdb2mat(['mimic3wdb/matched/', curr_rec_subfolder, '/', curr_rec(1:7), '/', curr_file], 1, 1, no_samps_to_download); % last sample, start at first sample
            else
                wfdb2mat(['mimic3wdb/matched/', curr_rec_subfolder, '/', curr_rec(1:7), '/', curr_file]); % last sample, start at first sample
            end
        end
        no_samps_downloaded = no_samps_downloaded + no_samps_to_download;

        % - extract information on this file
        [tm,signal,Fs,labels]=rdmat(filename);
        sig_names = extractfield(labels, 'Description');

        % - extract required data from this file
        if file_no == 1 || (subj_no == 35 && file_no == 2) || (subj_no == 22 && file_no == 5)

            % Add in data from the first file for this recording

            % - PPG
            rel_col = strcmp(sig_names, 'PLETH');
            data(inc_subj_no).ppg.v = signal(:,rel_col);
            data(inc_subj_no).ppg.fs = Fs;
            data(inc_subj_no).ppg.method = 'Fingertip PPG recorded using bedside monitor'; % From description at: https://physionet.org/content/mimicdb/1.0.0/

            % - ECG
            ecg_labels = {'II', 'I', 'V', 'III', 'MCL1'};
            % (manual intervention as the lead changed between files in these stays)
            if subj_no == 19 || subj_no == 35
                ecg_labels = {'V', 'II', 'I', 'III', 'MCL1'};
            end
            for ecg_label_no = 1 : length(ecg_labels)
                curr_ecg_label = ecg_labels{ecg_label_no};
                rel_col = strcmp(sig_names, curr_ecg_label);
                if sum(rel_col)
                    data(inc_subj_no).ekg.v = signal(:,rel_col);
                    data(inc_subj_no).ekg.fs = Fs;
                    data(inc_subj_no).ekg.method = ['ECG recorded using bedside monitor, lead ', curr_ecg_label];
                    break
                end
            end

            % - RESP
            rel_col = strcmp(sig_names, 'RESP');
            if ~isempty(rel_col)
                data(inc_subj_no).imp.v = signal(:,rel_col);
                data(inc_subj_no).imp.fs = Fs;
                data(inc_subj_no).imp.method = 'Impedance pneumography respiratory signal recorded at the chest using bedside monitor';
            end

            % - ABP
            rel_col = strcmp(sig_names, 'ABP');
            if ~isempty(rel_col)
                data(inc_subj_no).abp.v = signal(:,rel_col);
                data(inc_subj_no).abp.fs = Fs;
                data(inc_subj_no).abp.method = 'Invasive arterial blood pressure recorded using bedside monitor';
            end

            % - fixed data
            data(inc_subj_no).fix.subj_id = curr_subj;
            data(inc_subj_no).fix.rec_id = curr_rec;
            data(inc_subj_no).fix.files = possible_files;
            data(inc_subj_no).fix.af_status = stay_list.af(subj_no);

        else
            % Add in further data from additional files for this recording

            % - PPG
            rel_col = strcmp(sig_names, 'PLETH');
            data(inc_subj_no).ppg.v = [data(inc_subj_no).ppg.v; signal(:,rel_col)];

            % - ECG
            rel_col = strcmp(sig_names, curr_ecg_label);
            if ~sum(rel_col)
                error('Couldn''t find required ECG signal')
            end
            data(inc_subj_no).ekg.v = [data(inc_subj_no).ekg.v; signal(:,rel_col)];

            % - RESP
            rel_col = strcmp(sig_names, 'RESP');
            if ~isempty(rel_col)
                data(inc_subj_no).imp.v = [data(inc_subj_no).imp.v; signal(:,rel_col)];
                data(inc_subj_no).imp.fs = Fs;
                data(inc_subj_no).imp.method = 'Impedance pneumography respiratory signal recorded at the chest using bedside monitor';
            end

            % - ABP
            rel_col = strcmp(sig_names, 'ABP');
            if ~isempty(rel_col)
                data(inc_subj_no).abp.v = [data(inc_subj_no).abp.v; signal(:,rel_col)];
                data(inc_subj_no).abp.fs = Fs;
                data(inc_subj_no).abp.method = 'Invasive arterial blood pressure recorded using bedside monitor';
            end

        end

    end

    % check whether extra signals were present throughout the recording
    if length(data(inc_subj_no).imp.v) < length(data(inc_subj_no).ekg.v)
        data(inc_subj_no).imp = [];
    end
    if length(data(inc_subj_no).abp.v) < length(data(inc_subj_no).ekg.v)
        data(inc_subj_no).abp = [];
    end

    fprintf('\n Exported data lasting %d mins', round(length(data(end).ppg.v)/(60*data(end).ppg.fs)) );

    % - reduce duration to required duration
    durn = up.settings.req_durn;
    win_step = 1*60;
    clear s
    if (length(data(inc_subj_no).ppg.v)/data(inc_subj_no).ppg.fs) > durn
        
        % - identify period of data in which the PPG signal is most complete (and for which neither the PPG nor ECG contain nans)
        S = data(inc_subj_no).ppg;
        S = identify_periods_of_no_signal(S);
        S2 = data(inc_subj_no).ekg;
        S2 = identify_periods_of_no_signal(S2);
        win_starts_ppg = 1:(win_step*data(inc_subj_no).ppg.fs):(length(S.v)-(durn*data(inc_subj_no).ppg.fs));
        win_ends_ppg = win_starts_ppg + data(inc_subj_no).ppg.fs*durn;
        n_no_signal = nan(length(win_starts_ppg),1);
        for win_no = 1 : length(win_starts_ppg)
            % - identify window
            win_els = win_starts_ppg(win_no):win_ends_ppg(win_no);
            % - assess completeness of PPG data
            n_no_signal(win_no) = sum(S.no_signal(win_els)) + sum(S2.no_signal(win_els));
            % don't allow nans in either ECG or PPG or IMP - this step replaces the result with the max possible if there is at least one nan
            if sum(isnan(S.v(win_els))) || sum(isnan(S2.v(win_els))) || (~isempty(data(inc_subj_no).imp) && sum(isnan(data(inc_subj_no).imp.v(win_els))))
                n_no_signal(win_no) = n_no_signal(win_no) + 2*length(win_els);
            end
        end
        [~, most_complete_win] = min(n_no_signal);
        clear n_no_signal

        % - check that all sampling frequencies are the same
        fields = fieldnames(data);
        fs = [];
        for field_no = 1 : length(fields)
            eval(['curr_field_data = data(inc_subj_no).' fields{field_no} ';'])
            if ~isempty(curr_field_data) && sum(strcmp(fieldnames(curr_field_data), 'fs'))
                fs(end+1) = curr_field_data.fs;
            end
        end
        if length(unique(fs))>1
            error('The sampling frequencies of different waveforms should be the same')
        end

        % - specify which samples of signals to retain
        rel_win_start = win_starts_ppg(most_complete_win);
        rel_win_end = win_ends_ppg(most_complete_win);
        req_samps = rel_win_start:rel_win_end;
        clear win_starts* rel_win* S
    end
    clear win_step durn
    
    sigs = {'ekg', 'ppg', 'imp', 'abp'};
    for sig_no = 1 : length(sigs)
        eval(['rel_data = data(inc_subj_no).' sigs{sig_no} ';']);
        if ~isempty(rel_data)
            eval(['data(inc_subj_no).' sigs{sig_no} '.v = data(inc_subj_no).' sigs{sig_no} '.v(req_samps);']);
        end
    end
    clear req_samps

    fprintf('\n Refined to %d mins', round(length(data(end).ppg.v)/(60*data(end).ppg.fs)) );
end

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
    '\nMoody, B., Moody, G., Villarroel, M., Clifford, G. D., & Silva, I. (2020). MIMIC-III Waveform Database (version 1.0). PhysioNet. https://doi.org/10.13026/c2607m.' ...
    '\n\nThe MIMIC III Waveform Database is licensed under the ODbL 1.0 license.' ...
    '\n\nThe MIMIC-III database is described in the following publication:' ...
    '\nJohnson, A. E. W., Pollard, T. J., Shen, L., Lehman, L. H., Feng, M., Ghassemi, M., Moody, B., Szolovits, P., Celi, L. A., & Mark, R. G. (2016). MIMIC-III, a freely accessible critical care database. Scientific Data, 3, 160035. https://doi.org/10.1038/sdata.2016.35' ...
    '\n\nIt is available on PhysioNet, https://physionet.org/ :' ...
    '\nGoldberger, A., Amaral, L., Glass, L., Hausdorff, J., Ivanov, P. C., Mark, R., ... & Stanley, H. E. (2000). PhysioBank, PhysioToolkit, and PhysioNet: Components of a new research resource for complex physiologic signals. Circulation [Online]. 101 (23), pp. e215–e220.' ...
    '\n\n\nThe following annotations of AF and non-AF were used to create the dataset:', ...
    '\n\nBashar, Syed Khairul (2020): Atrial Fibrillation annotations of electrocardiogram from MIMIC III matched subset. figshare. Dataset. https://doi.org/10.6084/m9.figshare.12149091.v1', ...
    '\n\nBashar, S.K., Ding, E., Walkey, A.J., McManus, D.D. and Chon, K.H., 2019. Noise Detection in Electrocardiogram Signals for Intensive Care Unit Patients. IEEE Access, 7, pp.88357-88368. ﻿https://doi.org/10.1109/ACCESS.2019.2926199', ...
    '\n\nThis annotation information is reproduced under the terms of the CC BY 4.0 licence: https://creativecommons.org/licenses/by/4.0/'];

% This function defines which MIMIC III records should be downloaded. It uses
% information extracted from the following Excel file:
%
%   Bashar, Syed Khairul (2020): Atrial Fibrillation annotations of electrocardiogram from MIMIC III matched subset. figshare. Dataset. https://doi.org/10.6084/m9.figshare.12149091.v1
%
% Further information on the selection of records for download was provided in:
%
%   Bashar, S.K., Ding, E., Walkey, A.J., McManus, D.D. and Chon, K.H., 2019. Noise Detection in Electrocardiogram Signals for Intensive Care Unit Patients. IEEE Access, 7, pp.88357-88368. ﻿https://doi.org/10.1109/ACCESS.2019.2926199
%
% The information is reproduced under the terms of the CC BY 4.0 licence: https://creativecommons.org/licenses/by/4.0/

% identify AF and non-AF subjects
af_status = zeros(length(data),1);
for s = 1 : length(data)
    af_status(s) = data(s).fix.af_status;
end

% save data for AF and non-AF subjects in turn
fprintf('\n   - Saving data: ')
orig_data = data; clear data
% - AF
rel_subjs = af_status == 1;
data = orig_data(rel_subjs);
fprintf('AF (%d subjs)', length(data));
save([up.paths.local.root_folder, 'mimic_af_data'], 'data', 'source', 'license')
% - AF
rel_subjs = af_status == 0;
data = orig_data(rel_subjs);
fprintf(', non-AF (%d subjs)', length(data));
save([up.paths.local.root_folder, 'mimic_non_af_data'], 'data', 'source', 'license')

%% Closing messages

fprintf('\n\n - NB: this dataset also contains additional variables which have not been imported');

fprintf('\n\n ~~~ DONE ~~~')

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
curr_rec_subfolder = curr_rec(1:3);
curr_rec_hea_path = [up.paths.web.root_folder, curr_rec_subfolder, '/', curr_rec(1:7), '/', curr_rec, '.hea'];
rec_data = webread(curr_rec_hea_path);
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
rec_file_info.deb = datenum([temp{5}, ' ', temp{6}], 'HH:MM:SS.FFF dd/mm/yyyy');
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

function rec_file_info = get_rec_files_in_period(rec_file_info, onset_time, offset_time)

%% Calculate start and end of period
% - start time of relevant period for this subject
period_deb = rec_file_info.deb;
if ~isempty(onset_time)
    period_deb = period_deb + datenum(onset_time, 'HH:MM:SS');
end
% - start time of relevant period for this subject
if ~isempty(offset_time)
    period_fin = period_deb + datenum(offset_time, 'HH:MM:SS');
else
    period_fin = inf;
end

%% Identify recordings which are completely within this period
rec_file_info.in_period = false(size(rec_file_info.rec_deb));
for s = 1 : length(rec_file_info.rec_deb)

    if rec_file_info.rec_deb(s) >= period_deb && rec_file_info.rec_fin(s) <= period_fin
        rec_file_info.in_period(s) = 1;
    end

end

end

function S = identify_periods_of_no_signal(S)
% identifies periods when there is no signal (either flat line or nans)

uParams.analysis.durn_flat_line = 0.2; % threshold duration in secs, above which a flat line is considered to indicate no signal (below this it might just be due to a temporarily steady PPG value).

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

function records = define_stays_for_download

% This function defines which MIMIC III records should be downloaded. It uses
% information extracted from the following Excel file:
%
%   Bashar, Syed Khairul (2020): Atrial Fibrillation annotations of electrocardiogram from MIMIC III matched subset. figshare. Dataset. https://doi.org/10.6084/m9.figshare.12149091.v1
%
% Further information on the selection of records for download was provided in:
%
%   Bashar, S.K., Ding, E., Walkey, A.J., McManus, D.D. and Chon, K.H., 2019. Noise Detection in Electrocardiogram Signals for Intensive Care Unit Patients. IEEE Access, 7, pp.88357-88368. ﻿https://doi.org/10.1109/ACCESS.2019.2926199
%
% The information is reproduced under the terms of the CC BY 4.0 licence: https://creativecommons.org/licenses/by/4.0/
%

%% initial record selection
% For each subject, information has been extracted for one of their files:
%   - WF_FILENAME_1 (and 'Total time_1' and 'Onset Time_1' and 'Offset Time_1') was used if this had a duration of > 5 hours
%   - otherwise WF_FILENAME_2 was used
records.af = [0;0;1;1;0;1;1;0;0;0;1;0;0;0;0;1;0;1;1;0;0;1;1;0;0;0;0;1;0;1;1;1;1;1;1;1;1;0;1;0;0;1;1;0;1];
records.subj_id = [608;776;946;4490;4829;75796;9526;10391;13072;13136;14079;15852;16684;17344;19608;22954;23824;25117;26377;26964;29512;43613;50089;50384;55204;58932;62160;63039;63628;68956;69339;75371;77729;87275;79998;81349;85866;87675;89565;89964;92289;92846;94847;97547;99674];
records.file = {'p000608-2167-03-09-11-54';'p000776-2184-04-30-15-16';'p000946-2120-05-14-08-08';'p004490-2151-01-07-12-36';'p004829-2103-08-30-21-52';'p075796-2198-07-25-23-40';'p009526-2113-11-17-02-12';'p010391-2183-12-25-10-15';'p013072-2194-01-22-16-13';'p013136-2133-11-09-16-58';'p014079-2182-09-24-13-41';'p015852-2148-05-03-18-39';'p016684-2188-01-29-00-06';'p017344-2169-07-17-17-32';'p019608-2125-02-05-04-57';'p022954-2136-02-29-17-52';'p023824-2182-11-27-14-22';'p025117-2202-03-15-20-28';'p026377-2111-11-17-16-46';'p026964-2147-01-11-18-03';'p029512-2188-02-27-18-10';'p043613-2185-01-18-23-52';'p050089-2157-08-23-16-37';'p050384-2195-01-30-02-21';'p055204-2132-06-30-09-34';'p058932-2120-10-13-23-15';'p062160-2153-10-03-14-49';'p063039-2157-03-29-13-35';'p063628-2176-07-02-20-38';'p068956-2107-04-21-16-05';'p069339-2133-12-09-21-14';'p075371-2119-08-22-00-53';'p077729-2120-08-31-01-03';'p087275-2108-08-29-12-53';'p079998-2101-10-21-21-31';'p081349-2120-02-11-06-35';'p085866-2178-03-20-17-11';'p087675-2104-12-05-03-53';'p089565-2174-05-12-00-07';'p089964-2154-05-21-14-53';'p092289-2183-03-17-23-12';'p092846-2129-12-21-13-12';'p094847-2112-02-12-19-56';'p097547-2125-10-21-23-43';'p099674-2105-06-13-00-07'};
records.onset_time =  {'', '', '00:00:00', '00:00:00', '', '00:00:00', '00:00:28', '', '', '', '00:00:00', '', '', '', '', '00:00:00', '', '00:01:30', '00:03:41', '', '', '00:00:49', '00:06:29', '', '', '', '', '00:00:00', '', '00:00:03', '00:02:47', '00:28:07', '00:00:27', '00:00:17', '00:01:11', '00:00:00', '00:01:30', '', '00:00:51', '', '', '00:00:00', '00:01:14', '', '00:00:07'};
records.offset_time = {'', '', '06:13:12', '08:54:53', '', '42:56:15', '05:29:00', '', '', '', '03:57:49', '', '', '', '', '22:54:21', '', '12:38:58', '29:59:44', '', '', '30:16:16', '44:01:00', '', '', '', '', '30:20:34', '', '18:50:22', '08:39:47', '14:06:19', '15:20:01', '12:14:49', '20:16:36', '05:02:36', '23:38:35', '', '13:00:29', '', '', '12:17:40', '22:52:35', '', '11:53:00'};

%% refined record selection
use_refined_selection = 0;
if use_refined_selection
    % The record used for each subject was refined if the initial choice didn't result in any data.
    records.af = [0;0;1;1;0;1;1;0;0;0;1;0;0;0;0;1;0;1;1;0;0;1;1;0;0;0;0;1;0;1;1;1;1;1;1;1;1;0;1;0;0;1;1;0;1];
    records.subj_id = [608;776;946;4490;4829;75796;9526;10391;13072;13136;14079;15852;16684;17344;19608;22954;23824;25117;26377;26964;29512;43613;50089;50384;55204;58932;62160;63039;63628;68956;69339;75371;77729;87275;79998;81349;85866;87675;89565;89964;92289;92846;94847;97547;99674];
    records.file = {'p000608-2167-03-23-08-12';'p000776-2184-04-30-15-16';'p000946-2120-05-14-08-08';'p004490-2151-01-08-14-56';'p004829-2103-08-30-21-52';'p075796-2198-07-25-23-40';'p009526-2113-11-17-13-13';'p010391-2183-12-25-10-15';'p013072-2194-01-22-16-13';'p013136-2133-11-10-14-55';'p014079-2182-09-24-13-41';'p015852-2148-05-03-18-39';'p016684-2188-01-29-00-06';'p017344-2169-07-17-14-33';'p019608-2125-02-05-04-57';'p022954-2136-02-29-17-52';'p023824-2182-11-27-14-22';'p025117-2202-03-15-20-28';'p026377-2111-11-17-16-46';'p026964-2147-01-11-18-03';'p029512-2188-02-27-18-10';'p043613-2185-01-18-23-52';'p050089-2157-08-23-16-37';'p050384-2195-01-30-02-21';'p055204-2132-06-30-09-34';'p058932-2120-10-13-23-15';'p062160-2153-10-03-14-49';'p063039-2157-03-29-13-35';'p063628-2176-07-02-20-38';'p068956-2107-04-21-16-05';'p069339-2133-12-09-21-14';'p075371-2119-08-22-00-53';'p077729-2120-08-31-01-03';'p087275-2108-08-29-12-53';'p079998-2101-10-21-21-31';'p081349-2120-02-11-06-35';'p085866-2178-03-20-17-11';'p087675-2104-12-05-03-53';'p089565-2174-05-12-00-07';'p089964-2154-05-22-07-45';'p092289-2183-03-17-23-12';'p092846-2129-12-21-13-12';'p094847-2112-02-12-19-56';'p097547-2125-10-21-23-43';'p099674-2105-06-12-21-00'};
    records.onset_time =  {'', '', '00:00:00', '00:00:00', '', '00:00:00', '', '', '', '', '00:00:00', '', '', '', '', '00:00:00', '', '00:01:30', '00:03:41', '', '', '00:00:49', '00:06:29', '', '', '', '', '00:00:00', '', '00:00:03', '00:02:47', '00:28:07', '00:00:27', '00:00:17', '00:01:11', '00:00:00', '00:01:30', '', '00:00:51', '', '', '00:00:00', '00:01:14', '', ''};
    records.offset_time = {'', '', '06:13:12', '08:54:53', '', '42:56:15', '', '', '', '', '03:57:49', '', '', '', '', '22:54:21', '', '12:38:58', '29:59:44', '', '', '30:16:16', '44:01:00', '', '', '', '', '30:20:34', '', '18:50:22', '08:39:47', '14:06:19', '15:20:01', '12:14:49', '20:16:36', '05:02:36', '23:38:35', '', '13:00:29', '', '', '12:17:40', '22:52:35', '', ''};
    % - this didn't increase the number of included subjects, so I didn't use it
end

% Extract record ID
for s = 1 : length(records.file)
    records.stayID{s,1} = records.file{s}(1:7);
end

end
