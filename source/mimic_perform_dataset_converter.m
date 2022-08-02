function mimic_perform_dataset_converter
% MIMIC_PERFORM_DATASET_CONVERTER converts datasets into CSV and WFDB formats.
%   MIMIC_PERFORM_DATASET_CONVERTER converts MIMIC PERform datasets from
%   Matlab format into CSV and WFDB formats.
%
%   # Inputs
%   * specify the path of the dataset for conversion in the "up.matlab_dataset.path" variable within "setup_params" below.
%
%   # Outputs
%   * files : containing all the data written to new subfolders within the path specified by "up.matlab_dataset.path".
%
%   # Requirements:
%   The WFDB Toolbox is required to convert to WFDB format, which can be
%   downloaded from:
%   <https://physionet.org/physiotools/matlab/wfdb-app-matlab/>
%           
%   # Further Information
%   This version is provided to convert the MIMIC PERform datasets originally
%   described in:
%   Charlton P.H. et al. Detecting beats in the photoplethysmogram: benchmarking open-source algorithms, Physiological Measurement, 2022. <https://doi.org/10.1088/1361-6579/ac826d>
%
%   Further information on this study can be obtained at:
%   <https://peterhcharlton.github.io/project/ppg-beats/>
%
%   In addition, further information on the ppg-beats toolbox, including 
%   future versions, can be obtained at:
%   <https://github.com/peterhcharlton/ppg-beats>
%
%   # Feedback
%   Questions, comments, criticisms, complements, and contributions are welcome at:
%   <https://github.com/peterhcharlton/ppg-beats/issues>
%
%   # Version
%   v.1 - published on 1st August 2022 by Peter Charlton
%
%   # Source
%   Adapted from: RRest_dataset_converter.m, part of the RRest (<https://github.com/peterhcharlton/RRest/>) repository which is covered by the GNU public licence.
%   
%   # Author
%   Peter H. Charlton, University of Cambridge, June 2022.
%
%   # License - GPL-3.0
%      Copyright (c) 2022 Peter H. Charlton
%      This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%      This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%      You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

fprintf('\n\n~~~~~  Starting Dataset Conversion  ~~~~~');

%% setup parameters for waveform generation
up = setup_params;

%% Load data
[data, source, license] = load_data(up);

%% Convert data to WFDB format
if up.output.wfdb_format
    convert_to_wfdb(data, license, up);
end

%% Convert data to CSV format
if up.output.csv_format
    convert_to_csv(data, license, up);
end

end

function up = setup_params

fprintf('\n - Setting up Universal Parameters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% PARAMETERS TO BE SPECIFIED %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Path from which to load the dataset:
up.matlab_dataset.path = '/Users/petercharlton/Documents/Data/mimic_perform_af_dataset/mimic_perform_non_af_data.mat';
if ~exist(up.matlab_dataset.path, 'file'), warning('Specified dataset file does not exist. Please adjust ''up.matlab_dataset.path'' in ''setup_params'''); end
fprintf('\n    - Loading data from:\n        %s', up.matlab_dataset.path)

% folder containing the dataset
[folder,filename,ext] = fileparts(up.matlab_dataset.path);
up.matlab_dataset.folder = [folder, filesep];
up.matlab_dataset.filename = filename;
up.matlab_dataset.name = strrep(filename, '_data', '');

% direction of slash used on this computer:
up.slash = filesep;

% Type(s) of output
up.output.wfdb_format = true;
up.output.csv_format = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% OTHER PARAMETERS %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

up.csv.new_line = '\n';

% Check that the WFDB Matlab Toolbox has been installed
if ~exist('getWfdbClass', 'file')
    warning('Couldn''t find the WFDB Matlab Toolbox. Please install as described at the top of the file.')
end

end

function [data, source, license] = load_data(up)

fprintf('\n - Loading dataset in Matlab format');

load(up.matlab_dataset.path)

if ~exist('source', 'var')
    source = '';
end

if ~exist('license', 'var')
    license = '';
end

end

function convert_to_wfdb(data, license, up)

fprintf('\n - Converting to WFDB format');

% Check that the WFDB toolbox is installed
check_WFDB_installation;

% Create folder to save data in
data_type = 'wfdb';
save_folder = make_folder_to_save_data_in(data_type, up.matlab_dataset.name, up);

% add licence file
create_licence_file(license, save_folder)

% cycle through each subject
no_subjs = length(data);
file_names = cell(no_subjs,1);
for subj_no = 1 : no_subjs
    
    % create data matrix for this subject
    [data_mat, sig_names, fs, units, var_names] = create_data_matrix(data, subj_no, 0);
    
    % create filename for this subject
    file_names{subj_no} = create_subj_file_name(subj_no, up.matlab_dataset.name);
    
    % Obtain name of original file
    original_filename = obtain_original_filename(data, subj_no);
    
    % Obtain this subject's group
    group_name = obtain_subj_group(data, subj_no);

    % Create string of fixed variables for this subject
    descrip = ['<Original Subject ID>: ', data(subj_no).fix.subj_id, ' <Original Recording ID>: ', data(subj_no).fix.rec_id, ' <Original File>: ', original_filename, ' <subject group>: ' group_name];
    
    % convert this subject's data into WFDB format
    mat2wfdb(data_mat, file_names{subj_no}, fs, [], units, descrip, [], sig_names);
    
    clear data_mat fs units descrip sig_names var_names data_mat units
    
end

% create a RECORDS file
file_name = 'RECORDS';
file_path = [save_folder, file_name];
fid = fopen(file_path, 'w');
formatSpec = '%s\r\n';
for line_no = 1 : length(file_names)
    fprintf(fid, formatSpec, file_names{line_no});
end
fclose(fid)

% create an ANNOTATORS file
do_annotators = 0;
if do_annotators
    file_name = 'ANNOTATORS';
    file_path = [save_folder, file_name];
    fid = fopen(file_path, 'w');
    fclose(fid)
end

% create a DBS file
do_dbs = 0;
if do_dbs
    file_name = 'DBS';
    file_path = [save_folder, file_name];
    fid = fopen(file_path, 'w');
    fprintf(fid, [dataset_name, '\t', 'Synthetic ECG and PPG, modulated by respiration']);
    fclose(fid)
end

end

function convert_to_csv(data, license, up)

fprintf('\n - Converting to CSV format');

% Create folder to save data in
data_type = 'csv';
save_folder = make_folder_to_save_data_in(data_type, up.matlab_dataset.name, up);

% add licence file
create_licence_file(license, save_folder)

% cycle through each subject
no_subjs = length(data);
file_names = cell(no_subjs,1);
fprintf('\n - Writing CSV file for each subject:');
for subj_no = 1 : no_subjs
    
    fprintf('\n   - subj. %.3d', subj_no);

    % create data matrix for this subject
    [data_mat, sig_names, fs, units, var_names] = create_data_matrix(data, subj_no, 1);

    % limit precision of the data
    limit_precision = 0;
    if limit_precision
        data_mat = limit_precision(data_mat, var_names);
    end

    % create filename for this subject
    file_names{subj_no} = create_subj_file_name(subj_no, up.matlab_dataset.name);

    % Obtain name of original file
    original_filename = obtain_original_filename(data, subj_no);
    
    % Obtain this subject's group
    group_name = obtain_subj_group(data, subj_no);
    
    % create list of signals
    sig_text = [];
    for sig_no = 1 : length(sig_names)
        sig_text = [sig_text, sig_names{sig_no} '; and '];
    end
    sig_text = sig_text(1:end-6);

    % create a file of fixed variables for this subject
    curr_filename = [save_folder, file_names{subj_no}, '_fix.txt'];
    fid = fopen(curr_filename, 'wt');
    text_to_write = [file_names{subj_no}, up.csv.new_line, ...
        'Signals: ' sig_text, up.csv.new_line, ...
        'Sampling frequency: ' num2str(fs) ' Hz', up.csv.new_line, ...
        'Original Subject ID: ' data(subj_no).fix.subj_id, up.csv.new_line, ...
        'Original Recording ID: ' data(subj_no).fix.rec_id, up.csv.new_line, ... 
        'Original file: ' original_filename, up.csv.new_line, ...
        'Group: ' group_name, up.csv.new_line, ...
        ];
    fprintf(fid, text_to_write); clear text_to_write
    
    % close file
    fclose(fid); clear fid
    
    % create a file of data for this subject
    curr_filename = [save_folder, file_names{subj_no}, '_data.csv'];
    T = array2table(data_mat);  % Adapted from: https://uk.mathworks.com/matlabcentral/answers/281150-writing-a-matrix-with-header-into-a-csv-file
    T.Properties.VariableNames = var_names;
    writetable(T, curr_filename); clear data_mat T
    
    clear fs units descrip sig_names var_names T curr_filename data_mat units
    
end

fprintf('\n - Converted to CSV format');

end

function create_licence_file(license, save_folder)

if ~isempty(license)
    license_text = license.details;
else
    license_text = standard_license_text;
end

licence_file_path = [save_folder, 'LICENSE'];
fid = fopen(licence_file_path, 'wt');
fprintf(fid, license_text);
fclose all;

end

function license_text = standard_license_text

license_text = ['This dataset is licensed under the Open Data Commons Open Database License v1.0 (ODbL 1.0 license).' ...
    '\nFurther details are available at: https://opendatacommons.org/licenses/odbl/summary/' ...
    '\n\nThis dataset is derived from the MIMIC III Waveform Database:' ...
    '\nMoody, B., Moody, G., Villarroel, M., Clifford, G. D., & Silva, I. (2020). MIMIC-III Waveform Database (version 1.0). PhysioNet. https://doi.org/10.13026/c2607m.' ...
    '\n\nThe MIMIC III Waveform Database is licensed under the ODbL 1.0 license.' ...
    '\n\nThe MIMIC-III database is described in the following publication:' ...
    '\nJohnson, A. E. W., Pollard, T. J., Shen, L., Lehman, L. H., Feng, M., Ghassemi, M., Moody, B., Szolovits, P., Celi, L. A., & Mark, R. G. (2016). MIMIC-III, a freely accessible critical care database. Scientific Data, 3, 160035. https://doi.org/10.1038/sdata.2016.35' ...
    '\n\nIt is available on PhysioNet, https://physionet.org/ :' ...
    '\nGoldberger, A., Amaral, L., Glass, L., Hausdorff, J., Ivanov, P. C., Mark, R., ... & Stanley, H. E. (2000). PhysioBank, PhysioToolkit, and PhysioNet: Components of a new research resource for complex physiologic signals. Circulation [Online]. 101 (23), pp. e215–e220.'];


end

function save_folder = make_folder_to_save_data_in(data_type, dataset_name, up)

save_folder = [up.matlab_dataset.folder, dataset_name, '_', data_type, up.slash];
if ~exist(save_folder, 'dir')
    fprintf('\n    - Making folder in which to save data at:\n    - %s', save_folder)
    mkdir(save_folder)
end
cd(save_folder)

end

function [data_mat, sig_names, fs, units, var_names] = create_data_matrix(data, subj_no, do_time_vector)

data_mat(:,1) = data(subj_no).ppg.v;
data_mat(:,2) = data(subj_no).ekg.v;
if ~isempty(data(subj_no).imp)
    data_mat(:,3) = data(subj_no).imp.v;
    do_resp = 1;
else
    do_resp = 0;
end

var_names = {'PPG','ECG'};
if do_resp
    var_names(end+1) = {'resp'};
end

if do_time_vector
    time_vector = [0:length(data(subj_no).ppg.v)-1]/data(subj_no).ppg.fs; time_vector = time_vector(:);
    data_mat = [time_vector, data_mat];
    var_names = ['Time', var_names];
end

sig_names = {['PPG, ' data(subj_no).ppg.method], ['ECG, ' data(subj_no).ekg.method]};
units = 'NU/mV';
if do_resp
    sig_names(end+1) = {['resp, ' data(subj_no).imp.method]};
    units = [units, '/pm'];
end

fs = data(subj_no).ekg.fs;

end

function file_name = create_subj_file_name(subj_no, matlab_dataset_name)

file_name = [matlab_dataset_name, '_', sprintf('%.3d', subj_no)];

end

function check_WFDB_installation

% Check whether the 'wrsamp' function is available
if ~exist('wrsamp.m', 'file')
    error('This code requires the Waveform Database Software Package (WFDB) for MATLAB and Octave, available from: https://doi.org/10.13026/6zcz-e163 . Please install it and check it is visible in the Matlab path.');
else
    fprintf('\n    - Detected that the Waveform Database Software Package (WFDB) for MATLAB and Octave is installed.')
end

end

function original_filename = obtain_original_filename(data, subj_no)

if sum(strcmp(fieldnames(data(subj_no).fix), 'file'))
    original_filename = data(subj_no).fix.file;
else
    original_filename = [];
    for file_no = 1 : length(data(subj_no).fix.files)
        original_filename = [original_filename, data(subj_no).fix.files{1} '; '];
    end
    original_filename = original_filename(1:end-2);
end

end

function group_name = obtain_subj_group(data, subj_no)

if sum(strcmp(fieldnames(data(subj_no).fix), 'group'))
    group_name = data(subj_no).fix.group;
else
    group_name = '';
end

end

function data_mat = limit_precision(data_mat, var_names)

for var_no = 1 : length(var_names)
    % extract data for this variable
    curr_var = data_mat(:,var_no);

    % decide whether or not to reduce it's precision (i.e. number of decimal places)
    curr_var_range = range(curr_var);
    if curr_var_range > 0.5
        reduce_precision = 1;
    else
        reduce_precision = 0;
    end
    
    % reduce precision
    if reduce_precision
        curr_var = round(curr_var*1000)/1000;
    end

    % store results
    data_mat(:,var_no) = curr_var;

end


end