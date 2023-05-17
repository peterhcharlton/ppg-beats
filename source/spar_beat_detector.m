function [peaks, onsets] = spar_beat_detector(sig,fs)
% SPAR_BEAT_DETECTOR  SPAR PPG beat detector.
%   SPAR_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
%   using the 'Symmetric Projection Attractor Reconstruction' beat detector
%   
%   # Inputs
%   
%   * sig : a vector of PPG values
%   * fs  : the sampling frequency of the PPG in Hz
%   
%   # Outputs
%   * peaks : indices of detected pulse peaks
%   * onsets : indices of detected pulse onsets
%   
%   # Reference
%   C. Pettit and P. Aston, 'Photoplethysmogram (PPG) Beat Detection Using Symmetric Projection Attractor Reconstruction,' [in preparation].
%   
%   # Author
%   * Callum Pettit and Philip Aston - wrote the code: 'SPARcycleTimesPPG' and 'find_acl'
%   * Peter H. Charlton - modified the code slightly (as indicated), but mainly just wrote this wrapper
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Version
%   1.1
%   
%   # Licence
%      SPAR PPG Beat Detector as [downloaded/received] is Copyright ©, 2021-2022 University of Surrey.
%       
%      All Rights Reserved - You are allowed to alter and use this code only for your own non-commercial research purposes and this excludes commercial use of any kind. You are not permitted to copy and distribute this code to anyone under any circumstances. All rights remain with the Copyright holders.
%       
%      If you wish to copy, distribute or make any commercial use of the Software please contact Surrey’s IP & Licensing team at techtransfer@surrey.ac.uk to request an appropriate licence.
%       
%      Please contact Philip Aston (P.Aston@surrey.ac.uk) with any technical questions or collaboration requests.

N = 7;
[beat_indices_all,cycle_times_all] = SPARcycleTimesPPG(sig,fs,N);
peaks = nan(length(beat_indices_all)-1,1);
for beat_no = 1 : length(beat_indices_all)-1
    curr_ind = beat_indices_all(beat_no);
    [~, temp] = max(sig(curr_ind:beat_indices_all(beat_no+1)));
    peaks(beat_no) = curr_ind+temp-1;
end

peaks = tidy_beats(peaks);

onsets = pulse_onsets_from_peaks(sig, peaks);

end

function [beat_indices_all,cycle_times_all] = SPARcycleTimesPPG(signal_all,fs,N)

% This function finds the beats and cycle times for (human) PPG data
%
% Inputs:
%   signal_all: PPG signal in vector form
%   fs:         Sampling frequency that the data was recorded at
%   N:          Embedding dimension for SPAR method
%
% Outputs:
%   beat_indices_all: indices of points nearest to w=0 on each cycle
%   cycle_times_all:  the times of each individual cycle
%
% Written by Callum Pettit and Philip Aston
% August 2021

% set parameters

%global acl
acl = [];

k = 1;                 % the projection - we only use the first one
pltACL = 0;            % 0 = no plot, 1 = plot ACL function
pltAttr = 0;           % 0 = no plot, 1 = plot the attractor
HRmin = 40;            % heart rate minimum
HRmax = 200;           % heart rate maximum             (edited by PC from the original value of 120)
window_length = 20;    % window length in seconds

% define a time vector and find the number of windows

signal_all_len = length(signal_all);
time_all = linspace(0,(signal_all_len-1)/fs,signal_all_len)';
num_windows = floor(time_all(end)/window_length);
if num_windows==0
    num_windows = 1;
end

% initialise variables

beat_indices_all = [];
cycle_times_all = [];
tau = 0;
last_intersection = 0;
acl_time = 0;

% go through the windows

for win=1:num_windows
    
    % define the window with some overlap on the left
    
    window_start = last_intersection-(N-1)*tau/fs-0.2*acl_time;
    if win~=num_windows
        window = window_start<=time_all & time_all<=window_length*win;
    else
        window = window_start<=time_all; % include the data to the end for the last window
    end
    
    time = time_all(window);
    signal = signal_all(window);
    sig_len = length(signal);
    location = find(window,1);
    
    % find the average cycle length and hence tau
    redo_acl = 0;   % Line added by PC
    if redo_acl || isempty(acl) || acl==0    % Edited by PC (by adding "redo_acl || isempty(acl) ||")
        acl = find_acl(signal,fs,HRmin,HRmax,pltACL,win,...
            window_length*(win-1),window_length*(win));
    end
    acl_time = acl/fs;     % acl in time units
    tau = round(acl/N);    % time delay tau=acl/N
    if tau<5
        warning('The time delay parameter is tau = %d datapoints',tau)
    end
    
    % find v_N,k and w_N,k
    % formulas from:
    % J.V. Lyle and P.J Aston, Symmetric Projection Attractor Reconstruction:
    % Embedding in Higher Dimensions
    
    xN = zeros(N,sig_len-(N-1)*tau);  % the delay coordinates
    for i=1:N
        xN(i,:) = signal(1+(N-i)*tau:sig_len-(i-1)*tau);
    end
    j = 0:N-1;
    coeffs = sqrt(2/N)*[-cos(2*pi*(j+1)*k/N); sin(2*pi*(j+1)*k/N)];
    vwNk = coeffs*xN;
    v = vwNk(1,:)';    % v as a column vector
    w = vwNk(2,:)';    % w as a column vector
    length_v = length(v);
    vwTime = time(1+(N-1)*tau:sig_len);
    
    % rotate v and w by the optimal angle
    
    angle = pi*(0.5-k/N);    % optimal angle of rotation of the attractor
    vwRotated = [cos(angle) sin(angle); -sin(angle) cos(angle)] * [v w]';
    vRot = (vwRotated(1,:))';
    wRot = (vwRotated(2,:))';
    
    % plot the attractor and the section used
    
    if pltAttr==1
        maxvw = max([quantile(v,0.99),quantile(w,0.99)])*1.2;
        figure
        plot(v,w)
        hold on
        plot([0 2*maxvw*cos(angle)],[0 2*maxvw*sin(angle)],'r-')
        xlabel('v')
        ylabel('w')
        title(sprintf('The attractor for N = %d',N))
        axis([-maxvw,maxvw,-maxvw,maxvw],'square')
    end
    
    % detect crossings in one direction of the line w=0, v>0 in the
    % rotated coordinates
    
    Before = false(1,length_v);
    After = false(1,length_v);
    beats = false(1,length_v);
    for i = 1:length_v-1
        if wRot(i)>=0 && wRot(i+1)<0
            vIntersect = vRot(i)-(vRot(i)-vRot(i+1))*wRot(i)/(wRot(i)-wRot(i+1));
            if vIntersect>0  % check if v>0 at the intersection point
                Before(i) = true;
                After(i+1) = true;
                if wRot(i)<-wRot(i+1)
                    beats(i) = true;
                else
                    beats(i+1) = true;
                end
            end
        end
    end
    
    %%%%%%% Added by PC %%%%%%%
    do_pc_addition = 0;
    if do_pc_addition
        if sum(Before) == 0
            Before = false(1,length_v);
            After = false(1,length_v);
            beats = false(1,length_v);
            for i = 1:length_v-1
                if wRot(i)>=0 && wRot(i+1)<0
                    vIntersect = vRot(i)-(vRot(i)-vRot(i+1))*wRot(i)/(wRot(i)-wRot(i+1));
                    if vIntersect<0  % check if v<0 at the intersection point %%%%%% This is the change by PC
                        Before(i) = true;
                        After(i+1) = true;
                        if wRot(i)<-wRot(i+1)
                            beats(i) = true;
                        else
                            beats(i+1) = true;
                        end
                    end
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    wBefore = wRot(Before);
    wAfter = wRot(After);
    timeBefore = vwTime(Before);
    timeAfter = vwTime(After);
    
    % find the intersection times with w=0 by linear interpolation and
    % hence the cycle times
    
    intersections = (timeBefore-timeAfter).*(-wAfter)./(wBefore-wAfter) + timeAfter;
    %%% Added by PC %%%
    if isempty(intersections)
        continue
    end
    %%%%%%%%%%%%%%%%%%%
    last_intersection = intersections(end);
    cycle_times = diff(intersections)';
    cycle_times_all = [cycle_times_all cycle_times];
    
    % find the beat indices and correct for any mismatch between windows
    
    beat_indices = find(beats)+(N-1)*tau+location-1;
    if win~=1
        beat_correct = beat_indices_all(end)-beat_indices(1);
        beat_indices = beat_indices+beat_correct;
    end
    if win==1
        beat_indices_all = [beat_indices_all beat_indices];
    else
        beat_indices_all = [beat_indices_all beat_indices(2:end)];
    end
    
end

% fill in missed beats with estimated values

median_CT = median(cycle_times_all,'omitnan');
bias = 0.35;
upper_threshold = (1+bias)*median_CT;
lower_threshold = 0.7*median_CT;
missed_beats = find(cycle_times_all>upper_threshold);
num_missed_beats = length(missed_beats);
if num_missed_beats>0
    for j=num_missed_beats:-1:1
        consec_missed_beats = floor((cycle_times_all(missed_beats(j))-bias*median_CT)/median_CT);
        if missed_beats(j)==1 || missed_beats(j)==length(cycle_times_all) ||...
                (j~=num_missed_beats && missed_beats(j)+1==missed_beats(j+1)) ||...
                (j~=1 && missed_beats(j)-1==missed_beats(j-1)) ||...
                cycle_times_all(missed_beats(j)-1)<lower_threshold ||...
                cycle_times_all(missed_beats(j)+1)<lower_threshold
            extracycles = repmat(cycle_times_all(missed_beats(j))...
                /(consec_missed_beats+1),1,consec_missed_beats+1);
            extrabeats = beat_indices_all(missed_beats(j))...
                +round((beat_indices_all(missed_beats(j)+1)...
                -beat_indices_all(missed_beats(j)))...
                /(consec_missed_beats+1)*(1:consec_missed_beats));
        else
            slope = (cycle_times_all(missed_beats(j)+1)-cycle_times_all(missed_beats(j)-1))...
                /(consec_missed_beats+2);
            extracycles = cycle_times_all(missed_beats(j)-1)+slope*(1:consec_missed_beats+1);
            extrabeats = zeros(1,consec_missed_beats);
            for k=1:consec_missed_beats
                extrabeats(k) = beat_indices_all(missed_beats(j))...
                    +round(sum(extracycles(1:k)*fs));
            end
        end
        if missed_beats(j)==1
            cycle_times_all = [extracycles cycle_times_all(missed_beats(j)+1:end)];
        elseif missed_beats(j)==length(cycle_times_all)
            cycle_times_all = [cycle_times_all(1:missed_beats(j)-1) extracycles];
        else
            cycle_times_all = [cycle_times_all(1:missed_beats(j)-1)...
                extracycles cycle_times_all(missed_beats(j)+1:end)];
        end
        beat_indices_all = [beat_indices_all(1:missed_beats(j))...
            extrabeats beat_indices_all(missed_beats(j)+1:end)];
    end
end

% combine extra beats

extra_beats = find(cycle_times_all<lower_threshold);
num_extra_beats=length(extra_beats);
if num_extra_beats>0
    for j=num_extra_beats:-1:1
        if extra_beats(j)~=0
            if extra_beats(j)==1
                sumcycles = cycle_times_all(extra_beats(j)+1)+cycle_times_all(extra_beats(j));
                delta = 1;
            elseif extra_beats(j)==length(cycle_times_all)
                sumcycles = cycle_times_all(extra_beats(j)-1)+cycle_times_all(extra_beats(j));
                delta = -1;
            elseif j~=1 && extra_beats(j-1)==extra_beats(j)-1   % two neighbouring short beats
                sumcycles = cycle_times_all(extra_beats(j)-1)+cycle_times_all(extra_beats(j));
                delta = -1;
                if sumcycles<(1+bias)*median_CT && cycle_times_all(extra_beats(j)-1)>=lower_threshold
                    extra_beats(j-1) = 0;
                end
            elseif cycle_times_all(extra_beats(j)-1)<cycle_times_all(extra_beats(j)+1)
                sumcycles = cycle_times_all(extra_beats(j)-1)+cycle_times_all(extra_beats(j));
                delta = -1;
            else
                sumcycles = cycle_times_all(extra_beats(j)+1)+cycle_times_all(extra_beats(j));
                delta = 1;
            end
            if sumcycles<(1+bias)*median_CT
                cycle_times_all(extra_beats(j)+delta) = sumcycles;
                cycle_times_all(extra_beats(j)) = [];
                beat_indices_all(extra_beats(j)) = [];
            end
        end
    end
end

%%%%%%%%% Added by PC %%%%%%%%
beat_indices_all = sort(beat_indices_all);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function acl = find_acl(data,fs,HRmin,HRmax,plt,win,t0,t1)

% This functions finds the average cycle length (acl) using autocorrelation
%
% Inputs:
%   data: a column vector of data
%   fs:   sampling frequency of the data
%   HRmin,HRmax: parameters that specify a range of heart rate values to search for a maximum
%   plt:  0 = no plot; 1 = plot of the autocorrelation function and the maximum found
%   win:  window number - only used in the error message
%   t0,t1: time interval of the window - only used in the error message
%
% Output:
%   acl:  average cycle length
%
% Philip Aston
% August 2021

Tmin = floor(60*fs/HRmax);   % minimum possible cycle length
Tmax = ceil(60*fs/HRmin);    % maximum possible cycle length
Tmaxmum  =round(Tmax*1.1);

% check for too many NaN's in the data
% require at least 25% valid data to continue

tot = sum(isnan(data));
if tot/length(data)>0.25
    error('Too many NaNs in the data: %d%\n',round(100*tot/length(data)))
end

% filter the data to remove baseline wander using the filter recommended in:
% Yongbo Liang, Mohamed Elgendi, Zhencheng Chen and Rabab Ward 
% An optimal filter for short photoplethysmogram signals
% Scientific Data volume 5, Article number: 180076 (2018)

len = length(data);
data_ext = wextend(1,'sp0',data,len,'b'); % extend the data at both ends with a constant value
high_pass = 0.5;   % choose this parameter to get rid of low frequencies
order = 4;
[A,B,C,D] = cheby2(order,20,2*high_pass/fs,'high');
[filter_SOS,g] = ss2sos(A,B,C,D);
clean_data = filtfilt(filter_SOS,g,data_ext);
clean_data = clean_data(len+1:2*len);

% plot the original and the filtered data

% figure
% plot(data)
% hold on
% plot(clean_data,'r')
% pause
% close

% find the autocorrelation function and its peaks

acf = autocorr(clean_data,'NumLags',Tmaxmum);
[pk,loc] = findpeaks(acf(Tmin:Tmax),Tmin:Tmax);

if ~isempty(pk)
    [acfmax,m] = max(pk);  % find the highest peak
    acl = loc(m)-1;
    % check for a double peak at approximately half the cycle length
    if length(loc)>1 && m>1
        for j=1:m-1
            if pk(j)>=0.8*pk(m) ...
                    && ((0.45*loc(m)<=loc(j) && loc(j)<=0.55*loc(m))...
                    || (0.3*loc(m)<=loc(j) && loc(j)<=0.367*loc(m)))
                acfmax = pk(j);
                acl = loc(j)-1;
                break
            end
        end
    end
end

% produce a plot to show the maximum found or if no peaks are found

if plt==1 || isempty(pk)
    figure
    plot([Tmin Tmin],[-1 1],'m--')
    hold on
    plot([Tmax Tmax],[-1 1],'m--')
    text(Tmin+Tmaxmum/100,-0.9,'T_{min}')
    text(Tmax-Tmaxmum/12,-0.9,'T_{max}')
    plot(0:Tmaxmum,acf,'b')
    if ~isempty(pk)
        plot(acl,acfmax,'o','MarkerSize',7,'MarkerFaceColor','g','MarkerEdgeColor','g');
        title({['T = ' num2str(acl) ' (data points), T = ' num2str(acl/fs,3) ' (time units),', ...
        ' HR = ' num2str(60*fs/acl,3) ' (bpm)']},'FontSize',9,'FontWeight','normal');
    end
    axis([1 Tmaxmum -1 1])
    xlabel('Time Shift (data points)')
    ylabel('Autocorrelation function')
    hold off
end

if isempty(pk)
    error('No average cycle length found\nWindow: %d, start time: %f, end time: %f',...
        win,t0,t1)
end

end