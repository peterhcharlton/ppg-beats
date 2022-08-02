function [peaks, onsets] = heartpy_beat_detector(sig,fs)
% HEARTPY_BEAT_DETECTOR  HEARTPY PPG beat detector.
%   HEARTPY_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
%   using the 'HeartPy' beat detector
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
%   P. van Gent et al., 'HeartPy: A novel heart rate algorithm for the analysis of noisy signals,' Transportation Research Part F: Traffic Psychology and Behaviour, vol. 66, pp. 368-378, 2019. https://doi.org/10.1016/j.trf.2019.09.015>
%   
%   # Author
%   Peter H. Charlton, University of Cambridge, February 2022.
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Version
%   1.0
%   
%   # Source
%   This is a Matlab implementation of the original Python algorithm, described in:
%   van Gent P. et al., Analysing noisy driver physiology real-time using off-the-shelf sensors: Heart rate analysis software from the taking the fast lane project. J. Open Res. Softw. 2019, 7, <https://doi.org/10.5334/jors.241>
%   It was created based on v.1.2.5 of HeartPy downloaded on 7-June-2021 from: <https://github.com/paulvangentcom/heartrate_analysis_python> 
%   NB: The original HeartPy algorithm contains several extra features (such as artifact detection), which are not included in this code.
%   
%   # MIT License
%      Original:
%      MIT License
%      Copyright (c) 2021 Paul van Gent
%      Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
%      The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
%      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%
%      The code has since been implemented in Matlab format, and modified, by Peter H. Charlton.

%% Setup
hrdata = sig;
sample_rate = fs;
iterations = 2;

%% Pre-processing
% Peak enhancement
hrdata = enhance_peaks(hrdata);

%% Peak detection
% Process data
working_data = process(hrdata, sample_rate);

%% Output peaks
% Only output those that are kept
peaks = working_data.peaklist(working_data.binary_peaklist);

peaks = tidy_beats(peaks);

onsets = pulse_onsets_from_peaks(sig, peaks);

end

function working_data = process(hrdata, sample_rate)

% adapted from the 'process' function within 'heartpy.py'

% constants
bpmmin=40;
bpmmax=180;
windowsize=0.75;
reject_segmentwise=false;

working_data.hr = hrdata;
working_data.sample_rate = sample_rate;

% Calculate moving average
rol_mean = rolling_mean(hrdata, windowsize, sample_rate);

% Identify peaks
working_data = fit_peaks(hrdata, rol_mean, sample_rate, bpmmin, bpmmax);

% Calculate peak-to-peak intervals
working_data = calc_rr(working_data.peaklist, sample_rate, working_data);

% Check peaks
working_data = check_peaks(working_data.RR_list, working_data.peaklist, working_data.ybeat, reject_segmentwise, working_data);

% reminaing functions appear to be analysing peak-to-peak intervals, rather than identifying peaks, so not implemented here.

end

function working_data = fit_peaks(hrdata, rol_mean, sample_rate, bpmmin, bpmmax)
% from the 'fit_peaks' function within 'peakdetection.py'

% moving average values to test
ma_perc_list = [5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 150, 200, 300];

[rrsd.rrsd, rrsd.bpm, rrsd.ma_perc] = deal([]);
[valid_ma.rrsd, valid_ma.ma_perc] = deal([]);

update_dict = true;

for ma_val_no = 1 : length(ma_perc_list)
    ma_perc = ma_perc_list(ma_val_no);

    working_data = detect_peaks(hrdata, rol_mean, ma_perc, sample_rate, update_dict);
    bpm = ((length(working_data.peaklist)/(length(hrdata)/sample_rate))*60);
    
    rrsd.rrsd(end+1,1) = working_data.rrsd;
    rrsd.bpm(end+1,1) = bpm;
    rrsd.ma_perc(end+1,1) = ma_perc;

end

for ma_val_no = 1 : length(ma_perc_list)
    if rrsd.rrsd(ma_val_no) > 0.1 && rrsd.bpm(ma_val_no) >= bpmmin && rrsd.bpm(ma_val_no) <= bpmmax
        valid_ma.rrsd(end+1,1) = rrsd.rrsd(ma_val_no);
	    valid_ma.ma_perc(end+1,1) = rrsd.ma_perc(ma_val_no);
    end
end

if ~isempty(valid_ma.rrsd)
    [~, min_rrsd_el] = min(valid_ma.rrsd);
    working_data.best = valid_ma.ma_perc(min_rrsd_el);
    working_data = detect_peaks(hrdata, rol_mean, working_data.best, sample_rate, update_dict, working_data);
else
    warning('\n----------------\nCould not determine best fit for given signal. Please check the source signal.\n Probable causes:\n- detected heart rate falls outside of bpmmin<->bpmmax constraints\n- no detectable heart rate present in signal\n- very noisy signal (consider filtering and scaling)\nIf you\''re sure the signal contains heart rate data, consider filtering and/or scaling first.\n----------------\n');
end

end

function working_data = check_peaks(rr_arr, peaklist, ybeat, reject_segmentwise, working_data)
% excludes outlying peaks
%  an implementation of 'check_peaks' from 'peakdetection.py' in HeartPy.

% define RR range as mean +/- 30%, with a minimum of 300
mean_rr = mean(rr_arr);
thirty_perc = 0.3 * mean_rr;
if thirty_perc <= 300
    upper_threshold = mean_rr + 300;
    lower_threshold = mean_rr - 300;
else
    upper_threshold = mean_rr + thirty_perc;
    lower_threshold = mean_rr - thirty_perc;
end

% identify peaks to exclude based on RR interval
rem_idx = find((rr_arr <= lower_threshold) | (rr_arr >= upper_threshold));

working_data.removed_beats = peaklist(rem_idx);
working_data.removed_beats_y = ybeat(rem_idx);
working_data.binary_peaklist = true(size(peaklist));
working_data.binary_peaklist(rem_idx) = false;

if reject_segmentwise
    working_data = check_binary_quality(peaklist, working_data.binary_peaklist, working_data);
end

working_data = update_rr(working_data);

end

function working_data = check_binary_quality(peaklist, binary_peaklist, working_data)

% Constant
max_rejects = 3;

% The whole segment is rejected as it contains more than the specified 3 rejections per 10 beats.
idx = 0;
working_data.rejected_segments = [];
for i = 1: floor(length(binary_peaklist) / 10)
    if sum(binary_peaklist(idx:idx + 10)==0) > maxrejects
        binary_peaklist(idx:idx + 10) = 0;
        if idx + 10 < length(peaklist)
            working_data.rejected_segments(end+1,:) = [peaklist(idx), peaklist(idx + 10)];
        else
            working_data.rejected_segments(end+1,:) = [peaklist(idx), peaklist(end)];
        end
    end
    idx = idx+10;
end

end

function working_data = update_rr(working_data)
% Function that updates RR differences and RR squared differences based on corrected RR list
%  an implementation of 'update_rr' from 'analysis.py' in HeartPy

% - NB: this hasn't been completely implemented - it doesn't contain 'mask' variables. (and needs checking)

rr_source = working_data.RR_list;
b_peaklist = working_data.binary_peaklist;

rel_rrs = (b_peaklist(1:end-1) + b_peaklist(2:end)) == 2;
rr_list = rr_source(rel_rrs);
rr_diff = diff(rr_list);
rr_sqdiff = rr_diff.^2;

working_data.RR_list_cor = rr_list;
working_data.RR_diff = rr_diff;
working_data.RR_sqdiff = rr_sqdiff;

end

function working_data = detect_peaks(hrdata, rol_mean, ma_perc, sample_rate, update_dict, working_data)

rmean = rol_mean;

rol_mean = rmean + ((rmean ./ 100) .* ma_perc);
mn = mean(rmean ./ 100) .* ma_perc;
rol_mean = rmean + mn;

peaksx = find(hrdata > rol_mean);
peaksy = hrdata(peaksx);

peakedges = [1; find(diff(peaksx) > 1); length(peaksx)];
peaklist = [];

for i = 1: length(peakedges)
    try
        y_values = peaksy(peakedges(i):peakedges(i+1));
        [~, max_el] = max(y_values);
        peaklist = [peaklist; peaksx( peakedges(i) + max_el - 1)];
    catch
        % do nothing
    end
end

working_data.peaklist = peaklist;

if update_dict
    
    working_data.ybeat = hrdata(peaklist);
    working_data.rolling_mean = rol_mean;
    working_data = calc_rr(working_data.peaklist, sample_rate, working_data);
    if ~isempty(working_data.RR_list)
        working_data.rrsd = std(working_data.RR_list);
    else
        working_data.rrsd = inf;
    end
    
end

end

function working_data = calc_rr(peaklist, sample_rate, working_data)
% calculate peak-to-peak intervals
%  an implementation of 'calc_rr' from 'analysis.py' in HeartPy

% delete first peak if within first 150ms (signal might start mid-beat after peak)
if length(peaklist) > 0
    if peaklist(1) <= ((sample_rate / 1000.0) * 150)
        peaklist(1) = [];
        working_data.peaklist = peaklist;
        working_data.ybeat(1) = [];
    end
end

rr_list = (diff(peaklist) / sample_rate) * 1000.0;
rr_indices = [peaklist(1:end-1), peaklist(2:end)];
rr_diff = diff(rr_list);
rr_sqdiff = rr_diff.^2;
working_data.RR_list = rr_list;
working_data.RR_indices = rr_indices;
working_data.RR_diff = rr_diff;
working_data.RR_sqdiff = rr_sqdiff;

end

function rol_mean = rolling_mean(hrdata, windowsize, sample_rate)
% from the 'rolling_mean' function within 'datautils.py'

% calculate rolling mean
rol_mean = mean(sliding_window(hrdata, round(windowsize*sample_rate)), 1);

% need to fill 1/2 windowsize gap at the start and end
n_missvals = round((length(hrdata) - length(rol_mean))/2);
missvals_a = ones(1, n_missvals)*rol_mean(1);
missvals_b = ones(1, n_missvals)*rol_mean(end);
rol_mean = [missvals_a, rol_mean, missvals_b];

rol_mean = rol_mean(:);

end

function window_mat = sliding_window(hrdata, win_samps)
% an implementation of a Python function

% assemble sliding window matrix
win_starts = win_samps:(length(hrdata)-win_samps+1);
window_mat = nan(win_samps,length(win_starts));
for s = 1 : length(win_starts)
    window_mat(:,s) = hrdata(win_starts(s):win_starts(s)+win_samps-1);
end

end

function hrdata = enhance_peaks(hrdata)
% Peak enhancement
% Implementation of the 'enhance_peaks' function from 'preprocessing.py'

% constants
iterations = 2;

hrdata = scale_data(hrdata);

for i = 1 : iterations
    hrdata = hrdata.^2;
    hrdata = scale_data(hrdata);
end

end

function data = scale_data(data)

% thresholds
lower = 0;
upper = 1024;

% calculate range
rng = max(data) - min(data);

% find minimum
minimum = min(data);

% normalise
data = (upper - lower) .* ((data - minimum) ./ rng) + lower;

end