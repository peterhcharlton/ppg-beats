function [peaks, onsets] = abd_beat_detector(sig, fs)
% ABD_BEAT_DETECTOR  ABD PPG beat detector.
%   ABD_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
%   using the 'Automatic Beat Detection' beat detector
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
%   Aboy M et al., An automatic beat detection algorithm for pressure signals. IEEE Trans Biomed Eng 2005; 52: 1662-70. <https://doi.org/10.1109/TBME.2005.855725>
%   
%   # Author
%   Peter H. Charlton: King's College London (August 2017), University of Cambridge (February 2022)
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Version
%   1.0
%   
%   # Source
%   This script contains items either copied or modified from the pulse-analyse
%   (<https://github.com/peterhcharlton/pulse-analyse>) and RRest 
%   (<http://github.com/peterhcharlton/RRest/>) repositories which are covered
%   by the GNU public licence.
%   
%   # License - GPL-3.0
%      Copyright (c) 2022 Peter H. Charlton
%      This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%      This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%      You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

% Changes from original algorithm description:
%  1) PSD calculation method may not be exactly the same
%  2) Not conducted on windows of 10 s
%  3) Band-pass filtering may not produce exactly the right cut-offs
%  4) Wasn't sure what the lower and upper HR values were, so used 30 and 200 bpm
%  5) Changed the proportion of the lower HR value at which to draw the lower cut-off (from 50% to 80%)
%  6) Changed the percentile threshold for identifying peaks in the derivative from 90% to 75%
%  7) Haven't implemented harmonic PSD
%  8) HR estimation only considers HRs within the range of plausible HRs

% inputs
x = sig;  % signal
up = setup_up_abd_algorithm; % settings
w = fs*10; % window length (number of samples)
win_starts = 1:round(0.8*w):length(x);
win_starts(win_starts>=length(x)-w+1) = [];
win_starts = [win_starts, length(x)+1-w];

% before pre-processing
px = DetectMaxima(x,0);  % detect all maxima
if isempty(px)
    peaks = [];
    return
end

% detect peaks in windows
all_p4 = [];
all_hr = nan(length(win_starts)-1,1);
for win_no = 1 : length(win_starts)-1
    curr_els = win_starts(win_no):win_starts(win_no)+w-1;
    curr_x = x(curr_els);
    y1 = Bandpass(curr_x, fs, 0.9*up.fl_hz, 3*up.fh_hz); % Filter no.1
    hr = EstimateHeartRate(y1, fs, up); % Estimate HR from weakly filtered signal
    all_hr(win_no) = hr;
    y2 = Bandpass(curr_x, fs, 0.9*up.fl_hz, 2.5*hr/60); % Filter no.2
    y2_deriv = EstimateDeriv(y2); % Estimate derivative from highly filtered signal
    p2 = DetectMaxima(y2_deriv,up.deriv_threshold); % Detect maxima in derivative
    % plot(x), hold on, plot(p2, x(p2), 'or')
    % plot(y2_deriv), hold on, thresh = prctile(y2_deriv, up.deriv_threshold); plot([0,length(y2_deriv)], thresh*[1,1])
    y3 = Bandpass(curr_x, fs, 0.9*up.fl_hz, 10*hr/60);
    p3 = DetectMaxima(y3,60); % Detect maxima in moderately filtered signal
    % plot(x), hold on, plot(p3, x(p3), 'or')
    p4 = find_pulse_peaks(p2,p3);
%     plot(curr_x), hold on, plot(p4, curr_x(p4), 'or')
    all_p4 = [all_p4;win_starts(win_no)+p4-1];
end

all_p4 = unique(all_p4);

if ~isempty(all_p4)
    [peaks, fn] = IBICorrect(all_p4, px, median(all_hr), fs, up);
    peaks = unique(peaks);
else
    peaks = all_p4;
end
% plot(x), hold on, plot(p, x(p), 'or')

% plot(y1), hold on, plot(p4, y1(p4), 'or')

peaks = tidy_beats(peaks);

onsets = pulse_onsets_from_peaks(sig, peaks);

end

function up = setup_up_abd_algorithm

% plausible HR limits
up.fl = 30; % lower bound for HR
up.fh = 200; % upper bound for HR
up.fl_hz = up.fl/60;
up.fh_hz = up.fh/60;

% Thresholds
up.deriv_threshold = 75;            % originally 90
up.upper_hr_thresh_prop = 2.25;     % originally 1.75
up.lower_hr_thresh_prop = 0.5;     % originally 0.75

% Other parameters
up.win_size = 10; % in secs

end

function mt = DetectMaxima(sig,percentile)

% Table VI pseudocode

tr = prctile(sig, percentile);
ld = length(sig);

m = 1+find(sig(3:end) < sig(2:end-1) & sig(1:end-2) < sig(2:end-1));
mt = m(sig(m)>tr);

end

function bpf_sig = Bandpass(sig, fs, lower_cutoff, upper_cutoff)

% Filter characteristics: Eliminate VLFs (below resp freqs): For 4bpm cutoff
up.paramSet.elim_vlf.Fpass = 1.3*lower_cutoff;  % in Hz
up.paramSet.elim_vlf.Fstop = 0.8*lower_cutoff;   % in Hz
up.paramSet.elim_vlf.Dpass = 0.05;
up.paramSet.elim_vlf.Dstop = 0.01;

% Filter characteristics: Eliminate VHFs (above frequency content of signals)
up.paramSet.elim_vhf.Fpass = 1.2*upper_cutoff;  % in HZ
up.paramSet.elim_vhf.Fstop = 0.8*upper_cutoff;  % in HZ
up.paramSet.elim_vhf.Dpass = 0.05;
up.paramSet.elim_vhf.Dstop = 0.03;

% perform BPF
s.v = sig;
s.fs = fs;
s_evlf = elim_vlfs_abd(s, up);
s_filt = elim_vhfs(s_evlf, up);
bpf_sig = s_filt.v;

end

function s_filt = elim_vlfs_abd(s, up)
%% Filter pre-processed signal to remove frequencies below resp
% Adapted from RRest

%% Eliminate nans
s.v(isnan(s.v)) = mean(s.v(~isnan(s.v)));

%% Make filter
flag  = 'scale';
[N,Wn,BETA,TYPE] = kaiserord([up.paramSet.elim_vlf.Fstop up.paramSet.elim_vlf.Fpass]/(s.fs/2), [1 0], [up.paramSet.elim_vlf.Dstop up.paramSet.elim_vlf.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% % Gives a -3 dB cutoff at ? Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 0.0266;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(fs/2);

try
    s_filt.v = filtfilt(AMfilter.numerator, 1, s.v);
    s_filt.v = s.v-s_filt.v;
catch
    s_filt.v = s.v;
end
s_filt.fs = s.fs;
end

function s_filt = elim_vhfs(s, up)
%% Filter signal to remove VHFs
% Adapted from RRest

s_filt.fs = s.fs;

%% Eliminate nans
s.v(isnan(s.v)) = mean(s.v(~isnan(s.v)));

%% Check to see if sampling freq is at least twice the freq of interest
if (up.paramSet.elim_vhf.Fpass/(s.fs/2)) >= 1
    % then the fs is too low to perform this filtering
    s_filt.v = s.v;
    return
end

%% Create filter
% parameters for the low-pass filter to be used
flag  = 'scale';
[N,Wn,BETA,TYPE] = kaiserord([up.paramSet.elim_vhf.Fstop up.paramSet.elim_vhf.Fpass]/(s.fs/2), [1 0], [up.paramSet.elim_vhf.Dstop up.paramSet.elim_vhf.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% % Gives a -3 dB cutoff at cutoff_freq Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 0.3355;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(s.fs/2);

%% Remove VHFs
s_dt=detrend(s.v);
s_filt.v = filtfilt(AMfilter.numerator, 1, s_dt);
end

function hr = EstimateHeartRate(sig, fs, up)

% Estimate PSD
blackman_window = blackman(length(sig), 'periodic');
[pxx,f] = periodogram(sig, blackman_window,length(sig), fs);
% [ph, fh] = harmonicPSD(pxx,f);
ph = pxx; fh = f;

% Extract HR
rel_els = fh >= up.fl_hz & fh <= up.fh_hz;
rel_p = ph(rel_els);
rel_f = fh(rel_els);
[~,max_el] = max(rel_p);
hr = rel_f(max_el)*60;

end

function deriv = EstimateDeriv(sig)

% Savitzky Golay
deriv_no = 1;
win_size = 5;
deriv = savitzky_golay_abd(sig, deriv_no, win_size);

end

function deriv = savitzky_golay_abd(sig, deriv_no, win_size)

%% assign coefficients
% From: https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter#Tables_of_selected_convolution_coefficients
% which are calculated from: A., Gorry (1990). "General least-squares smoothing and differentiation by the convolution (Savitzky?Golay) method". Analytical Chemistry. 62 (6): 570?3. doi:10.1021/ac00205a007.

switch deriv_no
    case 0
        % - smoothing
        switch win_size
            case 5
                coeffs = [-3, 12, 17, 12, -3];
                norm_factor = 35;
            case 7
                coeffs = [-2, 3, 6, 7, 6, 3, -2];
                norm_factor = 21;
            case 9
                coeffs = [-21, 14, 39, 54, 59, 54, 39, 14, -21];
                norm_factor = 231;
            otherwise
                error('Can''t do this window size')
        end
    case 1
        % - first derivative
        switch win_size
            case 5
                coeffs = -2:2;
                norm_factor = 10;
            case 7
                coeffs = -3:3;
                norm_factor = 28;
            case 9
                coeffs = -4:4;
                norm_factor = 60;
            otherwise
                error('Can''t do this window size')
        end
        
    case 2
        % - second derivative
        switch win_size
            case 5
                coeffs = [2,-1,-2,-1,2];
                norm_factor = 7;
            case 7
                coeffs = [5,0,-3,-4,-3,0,5];
                norm_factor = 42;
            case 9
                coeffs = [28,7,-8,-17,-20,-17,-8,7,28];
                norm_factor = 462;
            otherwise
                error('Can''t do this window size')
        end
        
    case 3
        % - third derivative
        switch win_size
            case 5
                coeffs = [-1,2,0,-2,1];
                norm_factor = 2;
            case 7
                coeffs = [-1,1,1,0,-1,-1,1];
                norm_factor = 6;
            case 9
                coeffs = [-14,7,13,9,0,-9,-13,-7,14];
                norm_factor = 198;
            otherwise
                error('Can''t do this window size')
        end
        
    case 4
        % - fourth derivative
        switch win_size
            case 7
                coeffs = [3,-7,1,6,1,-7,3];
                norm_factor = 11;
            case 9 
                coeffs = [14,-21,-11,9,18,9,-11,-21,14];
                norm_factor = 143;
            otherwise
                error('Can''t do this window size')
        end
        
    otherwise
        error('Can''t do this order of derivative')        
end

if rem(deriv_no, 2) == 1
    coeffs = -1*coeffs;
end

A = [1,0];
filtered_sig = filter(coeffs, A, sig);
s=length(sig);
half_win_size = floor(win_size*0.5);
deriv=[filtered_sig(win_size)*ones(half_win_size,1);filtered_sig(win_size:s);filtered_sig(s)*ones(half_win_size,1)];
deriv = deriv/norm_factor;

end

function p4 = find_pulse_peaks(p2,p3)
p4 = nan(length(p2),1);
for k = 1 : length(p2)
    rel_el = find(p3>p2(k),1);
    if ~isempty(rel_el)
        p4(k) = p3(rel_el);
    end
end
p4 = p4(~isnan(p4));
end

function [pc, fn] = IBICorrect(p, m, hr, fs, up)

% Correct peaks' location error due to pre-processing
pc = nan(length(p),1);
for k = 1 : length(p)
    [~,rel_el] = min(abs(m-p(k)));
    pc1(k,1) = m(rel_el);    
end

% Correct false positives
% - identify FPs
d = diff(pc1)/fs;  % interbeat intervals in secs
fp = find_reduced_IBIs(d, median(hr), up);
% - remove FPs
pc2 = pc1(~fp);

% Correct false negatives
d = diff(pc2)/fs;  % interbeat intervals in secs
fn = find_prolonged_IBIs(d, median(hr), up);

pc = pc1;

end

function fn = find_prolonged_IBIs(IBIs, med_hr, up)

IBI_thresh = up.upper_hr_thresh_prop*60/med_hr;
fn = IBIs > IBI_thresh;

end

function fp = find_reduced_IBIs(IBIs, med_hr, up)

IBI_thresh = up.lower_hr_thresh_prop*60/med_hr;
fp = IBIs < IBI_thresh;

end