function [agreed_beats_inds, qual] = detect_ecg_beats(ecg, fs, options, no_signal_vector)
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
%    - options.beat_detectors - The ECG beat detectors to be used. Provide a cell containing either 1 or 2 out of: {'gqrs', 'jqrs', 'rpeakdetect'}
%    - options.hpf            - A logical indicating whether (1) or not (0) to high-pass filter the ECG to remove baseline wander
%
%   * no_signal_vector : (optional) a vector of the same length as the ECG vector, which indicates whether (1) or not (0) there was a signal at each sample
%   
%   # Outputs
%   
%   * agreed_beats_inds : the indices of detected beats which both beat detectors agree on
%   * qual : a logical indicating whether (1) or not (0) beats agreed between beat detectors for each ECG sample 
%   * temporary file - this script creates a temporary file in the current directory
%   
%   # Exemplary usage:
%
%       options.win_durn = 20;
%       [agreed_beats_inds, qual] = detect_ecg_beats(ecg, fs, options)
%   
%   # Requirements:
%   
%   When using the gqrs beat detector, the Waveform Database Software Package (WFDB) for MATLAB and Octave is required, available from:
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
%   It also uses the 'rpeakdetect.m' script (under the GNU public licence):
%   - Link: <http://www.mit.edu/~gari/CODE/ECGtools/ecgBag/rpeakdetect.m>
%   - Author: Gari Clifford
%   
%   # Version
%   1.0
%   
%   # License - GPL-3.0
%      This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%      This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%      You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

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

% - check how many beat detectors have been specified
if length(options.beat_detectors)>2, error('This code isn''t able to use more than two beat detectors'); end

% - Display startup message
if options.start_up_message, display_startup_message(options); end

%% Remove baseline wander
if options.hpf
    ecg = hpf_ecg(ecg,fs);
end

%% Detect beats using each beat detector in turn

for beat_detector_no = 1 : length(options.beat_detectors)

    % - Identify this beat detector
    curr_beat_detector = options.beat_detectors{beat_detector_no};
    
    % - Detect beats using this beat detector
    beat_inds{beat_detector_no} = run_beat_detector(ecg, fs, curr_beat_detector, options);

end

%% Identify 'correct' beats

% - only do (most of) the remaining steps if more than one beat detector was used
if length(options.beat_detectors) > 1

    if options.verbose, fprintf('\n - Identifying correct beat detections: '), end

    for beat_detector_no = 1 : length(beat_inds)

        % - deem all beats to be of low quality if any of the beat detectors didn't detect any beats
        if isempty(beat_inds{beat_detector_no})
            beat_correct_log = false(size(beat_inds{1}));
        end

    end

    if ~exist("beat_correct_log", 'var')

        % - Calculate difference matrix using gqrs as the reference
        diff_matrix = repmat(beat_inds{2}, [1, length(beat_inds{1})]) - beat_inds{1}';
        % - Find minimum differences
        min_abs_diff = min(abs(diff_matrix), [], 2);
        % - Identify correctly identified beats
        beat_correct_log = min_abs_diff<(options.qrs_tol_window*fs);

    end

    % identify beats as those which were agreed on by both beat detectors
    agreed_beats_inds = beat_inds{2}(beat_correct_log);

else

    % - if there was only one beat detector
    agreed_beats_inds = beat_inds{1};

end

% - Remove any repeated beat detections
agreed_beats_inds = sort(agreed_beats_inds);
if ~isempty(agreed_beats_inds)
    repeated_beats = [0; diff(agreed_beats_inds)< round(options.qrs_tol_window*fs)];
    agreed_beats_inds = agreed_beats_inds(~repeated_beats);
    clear repeated_beats
end

% - only do (most of) the remaining steps if more than one beat detector was used
if length(options.beat_detectors) > 1

    if options.verbose, fprintf('%d out of %d beats detected by %s deemed to be correct', length(agreed_beats_inds), length(beat_inds{1}), options.beat_detectors{1}), end

    %% Assess quality of beat detections in windows

    if options.verbose, fprintf('\n - Assessing quality of beat detections in windows: '), end

    % - Identify window start times
    win_starts = 0:(options.win_durn-options.win_overlap):(((length(ecg)-1)/fs)-options.win_durn);
    % - Create time vector
    t = 0:(1/fs):((length(ecg)-1)/fs);
    % - Identify times of disagreement between beat detectors
    t_disagree_beats = t(beat_inds{2}(~beat_correct_log));
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
        curr_win_beat_inds = t(beat_inds{2}) >= curr_start & t(beat_inds{2}) <= curr_end;
        if sum(curr_win_disagree_inds) == 0 && sum(curr_win_beat_inds)
            qual(rel_els) = true;
            continue
        end

        clear curr_start curr_end rel_els curr_win_disagree_inds
    end

    if options.verbose, fprintf('%.1f%% of ECG signal deemed to be of high quality', 100*mean(qual)), end

end

% Output a dummy 'qual' variable if there was only 1 beat detector
if ~exist('qual', 'var')
    qual = [];
end


end

function ecg = hpf_ecg(ecg,fs)

%% Specify filter characteristics
filt_characteristics.Fpass = 1.02;  % in Hz
filt_characteristics.Fstop = 0.6;   % in Hz      % 0.6 - 1.02 gives -3dB cutoff of 40 bpm (i.e. 0.67 Hz)
filt_characteristics.Dpass = 0.05;
filt_characteristics.Dstop = 0.01;

%% Filter to remove low frequencies
S.fs = fs; S.v = ecg;
s_filt = elim_vlfs(S, filt_characteristics);

%% Output
ecg = s_filt.v;

end

function s_filt = elim_vlfs(s, filt_characteristics)
%% Filter pre-processed signal to remove frequencies below resp
% Adapted from RRest

%% Eliminate nans
s.v(isnan(s.v)) = mean(s.v(~isnan(s.v)));

%% Make filter
flag  = 'scale';
[N,Wn,BETA,TYPE] = kaiserord([filt_characteristics.Fstop filt_characteristics.Fpass]/(s.fs/2), [1 0], [filt_characteristics.Dstop filt_characteristics.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% % Gives a -3 dB cutoff at ? Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 4.435e-3;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(s.fs/2);

if length(s.v) > (length(AMfilter.numerator)-1)*3
    % - Tukey Window to avoid edge effects
    win_edge_durn = 0.5; % in secs
    prop_win = win_edge_durn*2/((length(s.v)-1)/s.fs);
    tw = tukeywin(length(s.v),prop_win); 
    s_filt.v = filtfilt(AMfilter.numerator, 1, s.v.*tw);
    s_filt.v = s.v-s_filt.v;
else
    s_filt.v = s.v;
end
s_filt.fs = s.fs;
end

function beat_inds = run_beat_detector(ecg, fs, curr_beat_detector, options)

%% Display message saying which beat detector is being used
if options.verbose, fprintf('\n - Detecting beats using the %s beat detector: ', curr_beat_detector), end
    
%% gqrs beat detector
if strcmpi(curr_beat_detector, 'gqrs')

    % - check WFDB toolbox installation
    check_WFDB_installation;

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
    beat_inds = rdann('temp', 'qrs');
    
    % - display message about temporary files
    if options.verbose
        fprintf('\n   - temporary files created at: %s%stemp', cd, filesep)
    end

end

%% rpeakdetect detector

if strcmpi(curr_beat_detector, 'rpeakdetect')

    % - segment data into windows (because sometimes the peak detector doesn't detect beats on a long segment of ECG data, but does on shorter sub-segments)
    win_durn = 20*fs;
    win_overlap = 5*fs;
    win_starts = 1:(win_durn-win_overlap):(length(ecg)-win_durn);
    win_starts(end+1) = length(ecg)-win_durn;
    win_starts = unique(win_starts);

    beat_inds = [];
    
    % - detect beats on each window
    for win_no = 1 : length(win_starts)

        % - identify signal elements in this window
        els = win_starts(win_no):win_starts(win_no)+win_durn;
        
        % - detect beats
        [~, ~, ~, win_beat_inds, ~, ~]  = rpeakdetect(ecg(els),fs);

        win_beat_inds = win_beat_inds + win_starts(win_no) - 1;

        beat_inds = [beat_inds, win_beat_inds];

    end
    beat_inds = beat_inds(:);

    % remove repeated beat detections (due to windowing)
    beat_inds = sort(beat_inds);
    if ~isempty(beat_inds)
        repeated_beats = [0; diff(beat_inds)< round(options.qrs_tol_window*fs)];
        beat_inds = beat_inds(~repeated_beats);
        clear repeated_beats
    end

end

%% JQRS detector

if strcmpi(curr_beat_detector, 'jqrs')

    % - Set parameters
    HRVparams.PeakDetect.REF_PERIOD = 0.250;
    HRVparams.PeakDetect.THRES = .6;
    HRVparams.Fs = fs;
    HRVparams.PeakDetect.fid_vec = [];
    HRVparams.PeakDetect.SIGN_FORCE = [];
    HRVparams.PeakDetect.debug = 0;
    HRVparams.PeakDetect.windows = 15;
    HRVparams.PeakDetect.ecgType = 'MECG';

    % - Detect beats
    beat_inds = run_qrsdet_by_seg(ecg,HRVparams);

end

%% End matter
beat_inds = beat_inds(:); % make into a column vector
if options.verbose, fprintf('%d beats detected', length(beat_inds)), end

end

function options = setup_options(options)

specified_options = fieldnames(options);
possible_options = {'win_durn', 'win_overlap', 'verbose', 'qrs_tol_window', 'start_up_message', 'beat_detectors', 'hpf'};

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
               default_val = 0.15;  % in secs
           case 'start_up_message'
               default_val = 1;    % i.e. display the start-up message
           case 'beat_detectors'
               default_val = {'gqrs', 'jqrs'}; % the ECG beat detectors to use
           case 'hpf'
               default_val = 0;
       end
       % store the default value for this option
        eval(['options.' curr_option ' = default_val;']);
    end
    
end

end

function check_WFDB_installation

% Check whether the 'wrsamp' function is available
if ~exist('wrsamp.m', 'file')
    error('This code requires the Waveform Database Software Package (WFDB) for MATLAB and Octave, available from: https://doi.org/10.13026/6zcz-e163 . Please install it and check it is visible in the Matlab path.');
else
    fprintf('\n    - Detected that the Waveform Database Software Package (WFDB) for MATLAB and Octave is installed.')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rpeakdetect.m                                                     %
% From: http://www.mit.edu/~gari/CODE/ECGtools/ecgBag/rpeakdetect.m %
% accessed on 13 June 2022                                          %
% reproduced under the GNU general public license                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [hrv, R_t, R_amp, R_index, S_t, S_amp]  = rpeakdetect(data,samp_freq,thresh,testmode)

% [hrv, R_t, R_amp, R_index, S_t, S_amp]  = rpeakdetect(data, samp_freq, thresh, testmode); 
% R_t == RR points in time, R_amp == amplitude
% of R peak in bpf data & S_amp == amplitude of 
% following minmum. sampling frequency (samp_freq = 256Hz 
% by default) only needed if no time vector is 
% specified (assumed to be 1st column or row). 
% The 'triggering' threshold 'thresh' for the peaks in the 'integrated'  
% waveform is 0.2 by default.  testmode = 0 (default) indicates
% no graphics diagnostics. Otherwise, you get to scan through each segment.
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

%%%%%%%%%%% make threshold default 0.2 -> this was 0.15 on MIT data 
if nargin < 4
   testmode = 0;
end
%%%%%%%%%%% make threshold default 0.2 -> this was 0.15 on MIT data 
if nargin < 3
   thresh = 0.2;
end
%%%%%%%%%%% make sample frequency default 256 Hz 
if nargin < 2
   samp_freq = 256;
   if(testmode==1)
       fprintf('Assuming sampling frequency of %iHz\n',samp_freq);
   end
end

%%%%%%%%%%% check format of data %%%%%%%%%%
[a b] = size(data);
if(a>b)
 len =a;
end
if(b>a)
 len =b;
end

%%%%%%%%%% if there's no time axis - make one 
if (a | b == 1);
% make time axis 
  tt = 1/samp_freq:1/samp_freq:ceil(len/samp_freq);
  t = tt(1:len);
  x = data;
end
%%%%%%%%%% check if data is in columns or rows
if (a == 2) 
  x=data(:,1);
  t=data(:,2); 
end
if (b == 2)
  t=data(:,1);
  x=data(:,2); 
end

%%%%%%%%% bandpass filter data - assume 256hz data %%%%%
 % remove mean
 x = x-mean(x);
 
 % FIR filtering stage
 bpf=x; %Initialise
if( (samp_freq==128) & (exist('filterECG128Hz')~=0) )
        bpf = filterECG128Hz(x); 
end
if( (samp_freq==256) & (exist('filterECG256Hz')~=0) )
        bpf = filterECG256Hz(x); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%% differentiate data %%%%%%%%%%%%%%%%%%%%%%%%%%%
 dff = diff(bpf);  % now it's one datum shorter than before

%%%%%%%%% square data    %%%%%%%%%%%%%%%%%%%%%%%%%%%
 sqr = dff.*dff;   %
 len = len-1; % how long is the new vector now? 

%%%%%%%%% integrate data over window 'd' %%%%%%%%%%%%%%%%%%%%%%%%%
 d=[1 1 1 1 1 1 1]; % window size - intialise
 if (samp_freq>=256) % adapt for higher sampling rates
   d = [ones(1,round(7*samp_freq/256))]; 
 end
 % integrate
 mdfint = medfilt1(filter(d,1,sqr),10);
 % remove filter delay for scanning back through ECG
 delay = ceil(length(d)/2);
 mdfint = mdfint(delay:length(mdfint));
%%%%%%%%% segment search area %%%%%%%%%%%%%%%%%%%%%%%
 %%%% first find the highest bumps in the data %%%%%% 
 max_h = max (mdfint(round(len/4):round(3*len/4)));

 %%%% then build an array of segments to look in %%%%%
 %thresh = 0.2;
 poss_reg = mdfint>(thresh*max_h);

%%%%%%%%% and find peaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%% find indices into boudaries of each segment %%%
 left  = find(diff([0 poss_reg'])==1); % remember to zero pad at start
 right = find(diff([poss_reg' 0])==-1); % remember to zero pad at end
 
 %%%% loop through all possibilities  
 for(i=1:length(left))
    [maxval(i) maxloc(i)] = max( bpf(left(i):right(i)) );
    [minval(i) minloc(i)] = min( bpf(left(i):right(i)) );
    maxloc(i) = maxloc(i)-1+left(i); % add offset of present location
    minloc(i) = minloc(i)-1+left(i); % add offset of present location
 end

 R_index = maxloc;
 R_t   = t(maxloc);
 R_amp = maxval;
 S_amp = minval;   %%%% Assuming the S-wave is the lowest
                   %%%% amp in the given window
 S_t   = t(minloc);

%%%%%%%%%% check for lead inversion %%%%%%%%%%%%%%%%%%%
 % i.e. do minima precede maxima?
 if (minloc(length(minloc))<maxloc(length(minloc))) 
  R_t   = t(minloc);
  R_amp = minval;
  S_t   = t(maxloc);
  S_amp = maxval;
 end

%%%%%%%%%%%%
hrv  = diff(R_t);
resp = R_amp-S_amp; 

%%%%%%%%%%%%%%%%%%%%
if (testmode~=0)
figure(1)
hold off
subplot(4,1,1)
plot(t,x);hold on;plot(t,bpf,'r')
title('raw ECG (blue) and zero-pahse FIR filtered ECG (red)')
ylabel('ECG')
hold off;
subplot(4,1,2)
plot(t(1:length(mdfint)),mdfint);hold on;
%plot(t(1:length(sqr)),sqr);hold on;
plot(t,max(mdfint)*bpf/(2*max(bpf)),'r')
plot(t(left),mdfint(left),'og')
plot(t(right),mdfint(right),'om')
title('Integrated data with scan boundaries over scaled ECG')
ylabel('Int ECG')
hold off;
subplot(4,1,3)
plot(t,bpf,'r');hold on;
plot(R_t,R_amp,'+k');
plot(S_t,S_amp,'+g');
title('ECG with R-peaks (black) and S-points (green) over ECG')
ylabel('ECG+R+S')
hold off;
subplot(4,1,4)
hold off
plot(R_t(1:length(hrv)),hrv,'r+')
hold on
title('RR intervals')
ylabel('RR (s)')
hold off
fprintf('Press any key for next block of data\n');
pause
end

end