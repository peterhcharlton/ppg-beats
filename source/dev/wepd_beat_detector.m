function [peaks, onsets] = wepd_beat_detector(sig,fs)
% WEPD_BEAT_DETECTOR  WEPD PPG beat detector.
%   WEPD_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
%   using the 'waveform envelope peak detection' (WEPD) beat detector
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
%   D. Han et al., 'A Real-Time PPG Peak Detection Method for Accurate Determination of Heart Rate during Sinus Rhythm and Cardiac Arrhythmia,' Biosensors, vol. 12, no. 2, p. 82, 2022. <https://doi.org/10.3390/bios12020082>
%   
%   # Author
%   * Peter H. Charlton
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Version
%   1.0
%
%   # MIT License
%      Copyright (c) 2022 Peter H. Charlton
%      Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
%      The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
%      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

beat_indices_all = detect_beats(sig, fs);

% locate peaks because the timings of onsets aren't always exact
onsets = nan(length(beat_indices_all)-1,1);
for beat_no = 1 : length(beat_indices_all)-1
    curr_ind = beat_indices_all(beat_no);
    [~, temp] = min(sig(curr_ind:beat_indices_all(beat_no+1)));
    onsets(beat_no) = curr_ind+temp-1;
end

onsets = tidy_beats(onsets);

peaks = pulse_peaks_from_onsets(sig, onsets);

end

function onsets = detect_beats(sig, fs)
% this function detects beats in sig using the WEPD method, as described in Han et al., 2022.

%% Setup
plausible_hrs = [30, 200];

%% Bandpass filter
% bandpass filter from 0.5 - 5 Hz using 5th order zero phase elliptic filter
ftype = 'bandpass';
n = 3;
Rp = 10;
Rs = 50;
Wp = [0.5, 5]/(fs/2);
[b,a] = ellip(n,Rp,Rs,Wp,ftype);
bpf = filtfilt(b,a,sig);
%freqz(b,a)

%% Moving average filtering and differentiation
% Moving average filter three times

% setup
a = bpf;
M = round(fs/10);
% first moving average
b = mov_avg(a, M);
M = round(fs/9);
% second moving average
b = mov_avg(b, M);
% calculate first-order difference of filtered signal
c = nan(size(sig));
c(1:end-1) = b(2:end) - b(1:end-1);
% third moving average
b = mov_avg(b, M);

%% Normalize signal
non_nans = ~isnan(c);
d = (c - mean(c(non_nans)))/std(c(non_nans));

%% Invert signal
% find minima and maxima
d_min = find_minima(d);
d_max = find_minima(-1*d);
% calculate HRs from each
d_min_hr = find_hr(d_min, fs);
d_max_hr = find_hr(d_max, fs);
% decide whether to use minima or maxima
[can_use_min, can_use_max] = deal(true);
% - check plausibility of HRs
if d_min_hr > plausible_hrs(2) || d_min_hr < plausible_hrs(1)
    can_use_min = false;
end
if d_max_hr > plausible_hrs(2) || d_max_hr < plausible_hrs(1)
    can_use_max = false;
end
% - see if one has significantly fewer beats
diff_hr = d_min_hr - d_max_hr;
if abs(diff_hr) > 10
    if d_min_hr < d_max_hr
        can_use_max = false;
    else
        can_use_min = false;
    end
end
% - choose the one with the sharper peaks
if can_use_max && can_use_min
    min_sharpness = mean([d(d_min-1)-d(d_min); d(d_min+1)-d(d_min)]);
    max_sharpness = mean([d(d_max)-d(d_max-1); d(d_min)-d(d_min+1)]);
    if min_sharpness > max_sharpness
        can_use_max = false;
    else
        can_use_min = false;
    end
end

%% Local minima
if can_use_max
    min_els = d_max;
    inv_log = true;
else
    min_els = d_min;
    inv_log = false;
end

% invert if required
if inv_log
    d = -1*d;
end

    
%% Use envelope to remove false positive beats

% cubic spline interpolation to find envelope
env1 = interp1(min_els, d(min_els), 1:length(d), 'spline');
env1 = env1(:);

% hilbert filter to find envelope
first_non_nan = find(~isnan(d),1);
d(1:first_non_nan-1) = d(first_non_nan);
last_non_nan = find(~isnan(d),1,'last');
d(last_non_nan+1:end) = d(last_non_nan);
fl = round(1.5*fs);
[~,env2] = envelope(d,fl,'analytic');

% find beats
res = range(d)/50; % this step is done because env2 doesn't quite line-up on the PPG signal
d_temp = d - rem(d,res);
env1 = env1 - rem(env1,res);
env2 = env2 - rem(env2,res);
beats = find(env1==env2 & env1==d_temp);

%% Eliminate overlapped beats
tol = 0.3; % in secs
tol_samps = round(tol*fs);
rel_beats = true(size(beats));
finished = false; beat_no = 1;
while ~finished
    overlapping_beats = abs(beats-beats(beat_no))<tol_samps;
    [~,rel_el] = min(d(overlapping_beats));
    rel_beats(overlapping_beats) = false;
    temp = find(overlapping_beats,rel_el); temp = temp(end);
    rel_beats(temp) = true;
    if sum(overlapping_beats)
        beat_no = find(overlapping_beats,1,'last')+1;
    else
        beat_no = beat_no+1;
    end
    if beat_no >= length(beats)
        finished = true;
    end
end

rel_beats = beats(rel_beats);

%% Output beat locations
if inv_log
    peaks = rel_beats;
    onsets = pulse_onsets_from_peaks(d, peaks);
else
    onsets = rel_beats;
end

end

function b = mov_avg(a, M)
% calculate moving average

% setup
N = length(a);
b = nan(size(a));

% calculate for each index
for i = M : (N-M-1)
    % skip if this can't be calculated
    if i-M<1 || i+M>N
        continue
    end
    % calculate moving average for this index
    b(i-M+1) = (1/(2*M+1))*sum(a(i-M:i+M));
end
end

function min_els = find_minima(sig)
% find indices of local minima in a signal

min_els = 1+find( sig(2:end-1)<sig(1:end-2) & sig(2:end-1)<sig(3:end) );

end

function hr = find_hr(extrema_inds, fs)
% calculate heart rate from the indices of extrema (either maxima or minima)

hr = 60*length(extrema_inds)/(range(extrema_inds)/fs);

end