function collate_mimic_perform_ethnicity_dataset
% COLLATE_MIMIC_PERFORM_ETHNICITY_DATASET  Collates ethnicity dataset.
%   COLLATE_MIMIC_PERFORM_ETHNICITY_DATASET  Downloads and collates data from the MIMIC-III
%   Waveform Database Matched Subset from Black and White subjects, for PPG beat
%   detector performance evaluation.
%
%   # Inputs
%   * none (although see 'setup_up' for the parameters to be set)
%
%   # Working Files
%   The script downloads several files from PhysioNet, which are saved locally.
%
%   # Outputs
%   * files : Two MATLAB files containing the required data for PPG beat detector evaluation: 'mimic_B_data.mat', and 'mimic_W_data.mat'.
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

fprintf('\n ~~~ Downloading and collating MIMIC III \n     matched waveform database excerpt \n     containing data from Black and White patients. ~~~\n')

% - setup universal parameters
up = setup_up;

% - identify subjects for inclusion in the dataset
identify_subjs_for_dataset(up);

% - define stays for download
define_stays_for_download(up);

% - download dataset
download_dataset(up);

end

function up = setup_up

fprintf('\n - Setting up parameters')

% The following root folder contains subfolders for each subject (name
% 'S#' where # is the subject number)
up.paths.local.root_folder = '/Users/petercharlton/Documents/Data/mimiciii_ppg_ethnicity_beat_detection/';
if ~exist(up.paths.local.root_folder, 'dir'), warning('Specified folder does not exist.\n Please adjust ''up.paths.local.root_folder'' in ''setup_up'''); end

% Location of files online
up.paths.web.root_folder = 'https://physionet.org/files/mimic3wdb-matched/1.0/';
up.paths.web.records_numerics_file = [up.paths.web.root_folder, 'RECORDS-numerics'];
up.paths.web.records_file = [up.paths.web.root_folder, 'RECORDS'];
up.paths.web.admissions_csv_zipped = 'https://physionet.org/content/mimiciii/1.4/ADMISSIONS.csv.gz';

% Filepaths at which to save
up.paths.local.temp_folder = [up.paths.local.root_folder, 'downloaded_files', filesep];
if ~exist(up.paths.local.temp_folder, 'dir')
    fprintf('\n - Making directory in which to store downloaded data');
    mkdir(up.paths.local.temp_folder);
end

% Paths
up.paths.records_file_orig = [up.paths.local.root_folder, 'RECORDS'];
up.paths.admissions = [up.paths.local.root_folder, 'admissions.mat'];
up.paths.stay_list = [up.paths.local.root_folder, 'stay_list.mat'];
up.paths.admissions_csv_zipped = [up.paths.local.root_folder, 'ADMISSIONS.csv.gz'];
up.paths.admissions_file_orig = [up.paths.local.root_folder, 'ADMISSIONS.csv'];

% Settings
up.settings.req_durn = 1*60*10; % minimum duration of recording (in secs)
up.settings.no_subjs_per_ethnicity = 100;

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

function identify_subjs_for_dataset(up)

% download required files (downloads RECORDS file and admissions CSV)
download_req_files(up);

% identify matched records in MIMIC III (extracts subj_ids from RECORDS file)
subj_ids = extract_matched_recs(up);

% obtain ethnicities corresponding to these records (extracts ethnicities from admissions data and subj_ids)
eths = extract_ethnicities(subj_ids, up);

% find relevant subj IDs (either WHITE or BLACK ethnicity)
res = identify_rel_subj_ids(subj_ids, eths);

% extract ages of these subjects, and refine to only include adults (from 'patients' and icu_stays' data)
res = extract_age_and_genders(res, up);

% save lists of RECORDS corresponding to Black and White subjects
save_lists_of_records(res, up);

end

function download_req_files(up)

% Downloads the required files, and converts them to Matlab format:

fprintf('\n - Downloading required files from PhysioNet')

%% List of records
% - download text file containing list of records
if ~exist(up.paths.records_file_orig, 'file')
    outfilename = websave(up.paths.records_file_orig,up.paths.web.records_file);
end

%% Admissions data
if ~exist(up.paths.admissions, 'file')
    % - download zipped csv file containing admissions data
    if ~exist(up.paths.admissions_csv_zipped, 'file')
        msg_txt = ["You need to manually download the following file:"; up.paths.web.admissions_csv_zipped; "(NB: this is provided in the Command Window)"; "Then save it at:"; up.paths.admissions_csv_zipped; "Click 'OK' when this is done"];
        fprintf('\n --- The file can be downloaded from:\n     %s\n', up.paths.web.admissions_csv_zipped)
        f = msgbox(msg_txt,"Manual download required", "warn");
        waitfor(f);
    end
    if ~exist(up.paths.admissions_csv_zipped, 'file')
        error('The following file is required: %s\nIt can be downloaded from: %s', up.paths.admissions_csv_zipped, up.paths.web.admissions_csv_zipped)
    end
    % - unzip csv file
    gunzip(up.paths.admissions_csv_zipped,up.paths.local.root_folder);
    % - convert csv file to matlab format
    opts = delimitedTextImportOptions;
    opts.DataLines = [2, Inf];
    opts.VariableNames = ["ROW_ID", "SUBJECT_ID", "HADM_ID", "ADMITTIME", "DISCHTIME", "DEATHTIME", "ADMISSION_TYPE", "ADMISSION_LOCATION", "DISCHARGE_LOCATION", "INSURANCE", "LANGUAGE", "RELIGION", "MARITAL_STATUS", "ETHNICITY", "EDREGTIME", "EDOUTTIME", "DIAGNOSIS", "HOSPITAL_EXPIRE_FLAG", "HAS_CHARTEVENTS_DATA"];
    opts.VariableTypes = ["double", "double", "double", "datetime", "datetime", "datetime", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "datetime", "datetime", "categorical", "double", "double"];
    ADMISSIONS = readtable(up.paths.admissions_file_orig, opts);
    % - save in matlab format
    save(up.paths.admissions, 'ADMISSIONS');
end

end

function matched_subject_ids = extract_matched_recs(up)

% load data
fid = fopen(up.paths.records_file_orig);
matched_subject_ids = nan(100000,1);
counter = 0;
while ~feof(fid)
    counter = counter+1;
    curr_line = fgetl(fid);
    matched_subject_ids(counter) = str2double(curr_line(7:11));    % since of format 'p00/p000020/'
end
matched_subject_ids = matched_subject_ids(1:counter);
fclose(fid);

matched_subject_ids = sort(matched_subject_ids);

end

function eths = extract_ethnicities(subj_ids, up)

% load data
load(up.paths.admissions);

% find ethnicities corresponding to these subj ids
[~,admissions_rows,~] = intersect(ADMISSIONS.SUBJECT_ID, subj_ids);
eths = cellstr(ADMISSIONS.ETHNICITY(admissions_rows));

end

function res = identify_rel_subj_ids(subj_ids, eths)

white_els = contains(eths, 'WHITE');
black_els = contains(eths, 'BLACK');
els_to_include = white_els | black_els;

subj_id = subj_ids(els_to_include);
eth = eths(els_to_include);

res = table(subj_id, eth);

end

function res = extract_age_and_genders(res, up)

% only do this step if specified by the user (it requires additional files from PhysioNet)
refine_by_age = 0;
if ~refine_by_age
    return
end

% load data
load(up.paths.patients);
load(up.paths.icustays);

% extract gender of each patient
gender = cell(height(res),1);
for s = 1 : height(res)
    rel_el = PATIENTS.SUBJECT_ID == res.subj_id(s);
    gender(s) = cellstr(PATIENTS.GENDER(rel_el));
end
clear rel_el s
res.gender = gender; clear gender

% % Little experiment to see how many neonates there are (over 7,000)
% clear res
% subj_id = PATIENTS.SUBJECT_ID;
% res = table(subj_id);

% extract DOB of each patient
dob = nan(height(res),1);
for s = 1 : height(res)
    rel_el = PATIENTS.SUBJECT_ID == res.subj_id(s);
    dob(s) = datenum(PATIENTS.DOB(rel_el));
end
clear rel_el s

% extract INTIME of first icu stay for each patient - NB, this is helpful: https://physionet.org/content/mimic3wdb/1.0/RECORDS-neonates
intime = nan(height(res),1);
for s = 1 : height(res)
    rel_el = ICUSTAYS.SUBJECT_ID == res.subj_id(s);
    rel_el = find(rel_el,1);
    if ~isempty(rel_el)
        intime(s) = datenum(ICUSTAYS.INTIME(rel_el));
    end
end
clear rel_el s
res.age = (intime - dob)/365; % convert from days to years
clear intime dob

% eliminate subjects without an age
res = res(~isnan(res.age),:);

% eliminate subjects aged less than 18 years old
res = res(res.age>=18,:);

end

function save_lists_of_records(res, up)

% output numbers of subjects
white_els = contains(res.eth, 'WHITE');
black_els = contains(res.eth, 'BLACK');
% - totals
fprintf('\n Adult subjects in the matched subset with either Black or White ethnicity: %d', height(res));
fprintf('\n  - White: %d', sum(white_els));
fprintf('\n  - Black: %d', sum(black_els));

%% Save file with subject IDs of these subjects
for eth = {'Black', 'White'}
    curr_eth = eth{1,1};
    
    % - create string to write to text file
    eval(['rel_els = ' lower(curr_eth) '_els;']);
    rel_subj_id = res.subj_id(rel_els);
    str = '';
    for s = 1 : length(rel_subj_id)
        curr_record = num2str(rel_subj_id(s),'%06.f');
        curr_folder = ['p0', curr_record(2)];
        str = [str, curr_folder, '/', curr_record, '/', '\n'];
    end
    clear curr_record curr_folder rel_subj_id rel_els
    
    % - save to file
    filepath = [up.paths.local.root_folder, 'records_' curr_eth];
    fileid = fopen(filepath, 'w');
    fprintf(fileid, str);
    fclose(fileid);
    
end

end

function download_dataset(up)

%% Identify records which meet the requirements
% To be included, records must contain:
% - at least 1 hour continuously
% - of the required variables: (i) Pleth, (ii) ECG, and (iii) a respiration signal.

% load staylist
load(up.paths.stay_list)

%% Download data for each subject which meets the requirements

fprintf('\n - Downloading data for each subject which meets the requirements')
% cycle through each subject
inc_subj_no = 0;
eths = {};
for subj_no = 1 : length(stay_list.stayID)

    % check to see whether there are already enough of this subject's ethnicity
    curr_subj_eth = stay_list.eth{subj_no};
    no_of_this_eth = sum(strcmp(eths, curr_subj_eth));
    if no_of_this_eth >= up.settings.no_subjs_per_ethnicity
        continue
    end

    % extract current subject
    curr_subj = stay_list.stayID{subj_no};

    fprintf('\n   - Extracting data for subject no. %d (%s):', subj_no, curr_subj)

    % check whether this subject should be skipped (subjects which were manually identified as having unsuitable data)
    unsuitable_subjs = {'p000652', 'p014058', 'p012869'};
    if sum(strcmp(curr_subj, unsuitable_subjs))
        fprintf('\n      - unsuitable data')
        continue
    end

    % identify recording files for this subject
    curr_rec = stay_list.file{subj_no};
    rec_file_info = get_rec_file_info(curr_rec, up);
    if ~isstruct(rec_file_info)
        fprintf('\n      - couldn''t identify recording files in recording %s', curr_rec)
        continue
    end

    % identify recording files within the specified period
    rec_file_info = get_rec_files_in_period(rec_file_info);

    curr_rec_subfolder = curr_rec(1:3);

    % see whether each file for this recording contains the required signals
    possible_files = {}; curr_durn = 0;
    for file_no = 1 : min([25, length(rec_file_info.rec_deb)])

        fprintf(', file %d', file_no);

        curr_file = rec_file_info.rec_name{file_no};

        % - check whether this file is of a suitable duration
        file_durn = (rec_file_info.rec_no_samps(file_no)-1)/rec_file_info.fs; % in secs
        if file_durn < up.settings.req_durn
            fprintf(' (duration too short)');
            continue
        end
        if file_durn > 3*up.settings.req_durn
            continue
        end

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
        if ~sum(contains(vars, {'II', 'I', 'V', 'III', 'MCL1'})) || ~sum(contains(vars, 'PLETH')) || ~sum(contains(vars, 'RESP'))
            fprintf(' (didn''t have the signals)');
            possible_files = {};
            continue
        end

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
        fprintf('\n      - not enough data in recording %s', curr_rec)
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
    skipping_this_subj = false;
    for file_no = 1 : length(possible_files)

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
        try
            [tm,signal,Fs,labels]=rdmat(filename);
        catch
            fprintf(' (skipping this subject as couldn''t extract information from a file')
            inc_subj_no = inc_subj_no-1;
            data = data(1:inc_subj_no);
            skipping_this_subj = true;
            break
        end
        sig_names = extractfield(labels, 'Description');

        % - extract required data from this file
        if file_no == 1

            % Add in data from the first file for this recording

            % - PPG
            rel_col = find_rel_col(sig_names, 'PLETH');
            data(inc_subj_no).ppg.v = signal(:,rel_col);
            data(inc_subj_no).ppg.fs = Fs;
            data(inc_subj_no).ppg.method = 'PPG recorded using bedside monitor'; % From description at: https://physionet.org/content/mimicdb/1.0.0/

            % - ECG
            ecg_labels = {'II', 'I', 'V', 'III', 'MCL1'};
            for ecg_label_no = 1 : length(ecg_labels)
                curr_ecg_label = ecg_labels{ecg_label_no};
                rel_col = find_rel_col(sig_names, curr_ecg_label);
                if sum(rel_col)
                    data(inc_subj_no).ekg.v = signal(:,rel_col);
                    data(inc_subj_no).ekg.fs = Fs;
                    data(inc_subj_no).ekg.method = ['ECG recorded using bedside monitor, lead ', curr_ecg_label];
                    break
                end
            end

            % - RESP
            rel_col = find_rel_col(sig_names, 'RESP');
            if ~isempty(rel_col)
                data(inc_subj_no).imp.v = signal(:,rel_col);
                data(inc_subj_no).imp.fs = Fs;
                data(inc_subj_no).imp.method = 'Impedance pneumography respiratory signal recorded at the chest using bedside monitor';
            end

            % - ABP
            rel_col = find_rel_col(sig_names, 'ABP');
            if ~isempty(rel_col)
                data(inc_subj_no).abp.v = signal(:,rel_col);
                data(inc_subj_no).abp.fs = Fs;
                data(inc_subj_no).abp.method = 'Invasive arterial blood pressure recorded using bedside monitor';
            end

            % - fixed data
            data(inc_subj_no).fix.subj_id = curr_subj;
            data(inc_subj_no).fix.rec_id = curr_rec;
            data(inc_subj_no).fix.files = possible_files;
            data(inc_subj_no).fix.ethnicity = stay_list.eth{subj_no};

        else
            % Add in further data from additional files for this recording

            % - PPG
            rel_col = find_rel_col(sig_names, 'PLETH');
            data(inc_subj_no).ppg.v = [data(inc_subj_no).ppg.v; signal(:,rel_col)];

            % - ECG
            rel_col = find_rel_col(sig_names, curr_ecg_label);
            if ~sum(rel_col)
                error('Couldn''t find required ECG signal')
            end
            data(inc_subj_no).ekg.v = [data(inc_subj_no).ekg.v; signal(:,rel_col)];

            % - RESP
            rel_col = find_rel_col(sig_names, 'RESP');
            if ~isempty(rel_col)
                data(inc_subj_no).imp.v = [data(inc_subj_no).imp.v; signal(:,rel_col)];
                data(inc_subj_no).imp.fs = Fs;
                data(inc_subj_no).imp.method = 'Impedance pneumography respiratory signal recorded at the chest using bedside monitor';
            end

            % - ABP
            rel_col = find_rel_col(sig_names, 'ABP');
            if ~isempty(rel_col)
                data(inc_subj_no).abp.v = [data(inc_subj_no).abp.v; signal(:,rel_col)];
                data(inc_subj_no).abp.fs = Fs;
                data(inc_subj_no).abp.method = 'Invasive arterial blood pressure recorded using bedside monitor';
            end

        end

    end

    % check whether this subject is to be skipped
    if skipping_this_subj
        continue
    end

    % check whether extra signals were present throughout the recording
    if length(data(inc_subj_no).imp.v) < length(data(inc_subj_no).ekg.v)
        data(inc_subj_no).imp = [];
    end
    if length(data(inc_subj_no).abp.v) < length(data(inc_subj_no).ekg.v)
        data(inc_subj_no).abp = [];
    end

    fprintf('\n     - Exported data lasting %d mins', round(length(data(end).ppg.v)/(60*data(end).ppg.fs)) );


    % - reduce duration to 10 min if necessary
    durn = 10*60;
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

    fprintf('\n     - Refined to %d mins', round(length(data(end).ppg.v)/(60*data(end).ppg.fs)) );

    % Review this subject's data to see if it is suitable
    % - check for flat line in EKG or PPG
    waves_to_check = {'ekg', 'ppg'};
    skip_this_subj = false;
    for wave_no = 1 : length(waves_to_check)

        curr_wave = waves_to_check{wave_no};
        eval(['curr_data = data(inc_subj_no).' curr_wave '.v;']);

        % if there is a flat line for at least 10% of the signal, then remove this subject.
        if sum(curr_data == mode(curr_data)) > (0.1*length(curr_data))
            skip_this_subj = true;
            continue
        end

        % if there are any nans, then remove this subject
        if sum(isnan(curr_data))>0
            skip_this_subj = true;
            continue
        end
            
    end
    if skip_this_subj
        data(inc_subj_no) = [];
        inc_subj_no = inc_subj_no-1;
        fprintf('\n      - flat line in %s in recording %s', curr_wave, curr_rec)
        continue
    end

    % calculate number of subjects of each ethnicity
    eths = cell(length(data),1);
    for s = 1 : length(eths)
        eths{s} = data(s).fix.ethnicity;
    end
    fprintf('\n     - Current number of subjects in the dataset: %d', length(data))
    unique_eths = unique(eths);
    for s = 1 : length(unique_eths)
        curr_eth = unique_eths{s};
        fprintf('\n        - %s: %d', curr_eth, sum(strcmp(eths, curr_eth)));
    end

end
clear eths

%% Save data

% add source details
source.date_of_conversion = datetime('now');
source.matlab_conversion_script = [mfilename('fullpath'), '.m'];
source.raw_data_path = up.paths.local.temp_folder;

% identify subjects of different ethnicities
eths = cell(length(data),1);
for s = 1 : length(eths)
    eths{s} = data(s).fix.ethnicity;
end

% save data for subjects of different ethnicities in turn
fprintf('\n   - Saving data: ')
orig_data = data; clear data
% - cycle through ethnicities
unique_eths = unique(eths);
for s = 1 : length(unique_eths)
    curr_eth = unique_eths{s};
    fprintf([curr_eth, ', ']);
    rel_subjs = strcmp(eths, curr_eth);
    data = orig_data(rel_subjs);
    save([up.paths.local.root_folder, 'mimic_' curr_eth '_data'], 'data', 'source')
end

%% Closing messages

fprintf('\n\n - NB: this dataset also contains additional variables which have not been imported');

fprintf('\n\n ~~~ DONE ~~~')

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
curr_rec_subfolder = curr_rec(1:3);
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

% - skip if this has been done
if exist(up.paths.stay_list, 'file')
    return
end

fprintf('\n - Defining stays for download')

% This function defines which MIMIC III records should be downloaded. It uses
% information generated by the 'mimic_spo2_assessment' matlab script

% Extract relevant subject IDs from text file
fprintf('\n   - Extracting relevant subject IDs from the text file generated by ''mimic_spo2_assessment''');
[records.eth, records.path, records.subj_id] = deal(cell(0));
for eth = {'Black', 'White'}
    curr_eth = eth{1,1};
    curr_filepath = [up.paths.local.root_folder, 'records_' curr_eth];
    fid = fopen(curr_filepath);
    while ~feof(fid)
        records.path{end+1,1} = fgetl(fid);
        records.subj_id{end+1,1} = records.path{end}(5:end-1);
        records.eth{end+1,1} = curr_eth(1);
    end
    fclose(fid);
    clear fid curr_filepath curr_eth
end

% Identify corresponding records (just taking the first one for each subject)
fprintf('\n   - Identifying the first record for each subject')
% - download lsit of numerics records
numerics_records = webread(up.paths.web.records_numerics_file);
numerics_records = split(numerics_records, char(10)); % char(10) is the newline character
% - identify first numerics record for each subj
records.file = cell(size(records.subj_id));
for record_no = 1 : length(records.subj_id)
    curr_subj_id = records.subj_id{record_no};
    rel_numerics_record_el = find(contains(numerics_records, curr_subj_id),1);
    if ~isempty(rel_numerics_record_el)
        rel_numerics_file = numerics_records{rel_numerics_record_el};
        temp = strfind(rel_numerics_file, '/');
        records.file{record_no} = rel_numerics_file(temp(2)+1:end-1); % excluded the 'n' from the end
    end
end
% - exclude records without a file
els_to_exclude = cellfun(@isempty, records.file);
fields = fieldnames(records);
for field_no = 1 : length(fields)
    curr_field = fields{field_no};
    eval(['records.' curr_field ' = records.' curr_field, '(~els_to_exclude);']);
end

% Extract record ID
for s = 1 : length(records.file)
    records.stayID{s,1} = records.file{s}(1:7);
end

% save staylist
stay_list = records;
save(up.paths.stay_list, 'stay_list')

end
