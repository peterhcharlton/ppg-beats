function [peaks, onsets, mid_amps] = detect_ppg_beats(s, beat_detector)
% DETECT_PPG_BEATS  detects beats in PPG.
%   DETECT_PPG_BEATS detects beats in a photoplethysmogram (PPG) signal
%   using a specified beat detector
%   
%   # Inputs
%   
%   * s : a structure containing the following fields:
%    - v : a vector of PPG values
%    - fs : the sampling frequency of the PPG in Hz
%    
%   * beat_detector  - a string specifying the beat detector to be used
%   
%   # Outputs
%   * peaks : indices of pulse peaks
%   * onsets : indices of pulse onsets
%   * mid_amps : indices of mid-points on the upslope between onsets and peaks (defined as the mid-amplitude)
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Author
%   Peter H. Charlton, University of Cambridge, February 2022.
%   
%   # Version
%   1.1
%   
%   # License - GPL-3.0
%      Copyright (c) 2022 Peter H. Charlton
%      This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%      This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%      You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

%% Detect peaks and onsets in pulse signal

eval(['[peaks, onsets] = ' lower(beat_detector) '_beat_detector(s.v,s.fs);']);

%% Tidy up peaks and onsets

[peaks, onsets] = tidy_peaks_and_onsets(s.v, peaks,onsets);

%% Add mid-amplitude points
mid_amps = calc_mid_amp_points(s.v, peaks, onsets);

end

function [peaks, onsets] = tidy_peaks_and_onsets(sig, peaks,onsets)

% Tidy up peaks and onsets to make them conform to the following rules:
% (i) No two points at the same time
% (ii) At least one local minimum between consecutive peaks
% (iii) At least one local maximum between consecutive onsets
% (iv) Alternates between onsets and peaks
% (v) Starts with onset, and ends with peak

%% (i) No two points at the same time

% remove any repeated peaks (or onsets)
peaks = unique(peaks);
onsets = unique(onsets);

% If there is a peak and onset at the same index, then remove them both
repeated_vals = intersect(peaks,onsets);
peaks = setxor(peaks, repeated_vals);
onsets = setxor(onsets, repeated_vals);
clear repeated_vals

%% (ii) At least one local minimum between consecutive peaks

% If there are two peaks without a local minimum between them, then remove the lower one
local_min = islocalmin(sig);
finished = false;
while ~finished
    els_to_remove = [];
    for peak_no = 1 : length(peaks)-1
        rel_els = (peaks(peak_no)+1) : (peaks(peak_no+1)-1);
        if sum(local_min(rel_els)) == 0
            peak_vals = sig(peaks(peak_no:peak_no+1));
            [~, el_to_remove] = min(peak_vals);
            els_to_remove(end+1) = peak_no + el_to_remove - 1;
        end
    end 
    peaks(els_to_remove) = [];
    if isempty(els_to_remove)
        finished = true;
    end
end
clear local_min finished els_to_remove peak_no rel_els peak_vals el_to_remove

%% (iii) At least one local maximum between consecutive onsets

% If there are two onsets without a local maximum between them, then remove the higher one
local_max = islocalmax(sig);
finished = false;
while ~finished
    els_to_remove = [];
    for onset_no = 1 : length(onsets)-1
        rel_els = (onsets(onset_no)+1) : (onsets(onset_no+1)-1);
        if sum(local_max(rel_els)) == 0
            onset_vals = sig(onsets(onset_no:onset_no+1));
            [~, el_to_remove] = max(onset_vals);
            els_to_remove(end+1) = onset_no + el_to_remove - 1;
        end
    end 
    onsets(els_to_remove) = [];
    if isempty(els_to_remove)
        finished = true;
    end
end
clear local_max finished els_to_remove onset_no rel_els onset_vals el_to_remove

%% (iv) Alternates between onsets and peaks

% If there are two consecutive peaks, then insert an onset between them
peak_log = [true(size(peaks)); false(size(onsets))];
[els,order] = sort([peaks;onsets]);
peak_log = peak_log(order);
bad_els = find(diff(peak_log)==0 & peak_log(1:end-1)); % repeated peaks
if ~isempty(bad_els)  % if there is a repeated peak
    new_onsets = nan(length(bad_els),1);
    for bad_el_no = 1 : length(bad_els)   % cycle through each repeated peak
        curr_pks = [els(bad_els(bad_el_no)),els(bad_els(bad_el_no)+1)];
        [~,temp]=min(sig(curr_pks(1):curr_pks(2)));
        new_onsets(bad_el_no,1) = curr_pks(1) -1 + temp;
    end
    onsets = sort([onsets; new_onsets]);
end
clear new_onsets temp curr_pks bad_el_no bad_els peak_log order els

% If there are two consecutive onsets, then insert a peak between them
onset_log = [false(size(peaks)); true(size(onsets))];
[els,order] = sort([peaks;onsets]);
onset_log = onset_log(order);
bad_els = find(diff(onset_log)==0 & onset_log(1:end-1)); % repeated onsets
if ~isempty(bad_els)  % if there is a repeated peak
    new_peaks = nan(length(bad_els),1);
    for bad_el_no = 1 : length(bad_els)   % cycle through each repeated peak
        curr_trs = [els(bad_els(bad_el_no)),els(bad_els(bad_el_no)+1)];
        [~,temp]=max(sig(curr_trs(1):curr_trs(2)));
        new_peaks(bad_el_no,1) = curr_trs(1) -1 + temp;
    end
    peaks = sort([peaks; new_peaks]);
end
clear new_peaks temp curr_trs bad_el_no bad_els onset_log order els

%% (v) Starts with onset, and ends with peak

% Make sure that the first onset is before the first peak, and the last peak is after the last onset
finished = false;
while ~finished
    if ~isempty(onsets) && ~isempty(peaks) && onsets(1)>peaks(1)
        peaks(1) = [];
    else
        finished = true;
    end
end
finished = false;
while ~finished
    if ~isempty(peaks) && ~isempty(onsets) && peaks(end)<onsets(end)
        onsets(end) = [];
    else
        finished = true;
    end
end

% if no peaks (or onsets) were detected, then don't output any indices for either
if isempty(peaks)
    onsets = [];
end
if isempty(onsets)
    peaks = [];
end

end

function mid_amps = calc_mid_amp_points(sig, peaks, onsets)
% calculate mid-amplitude points between onsets and peaks

mid_amps = nan(length(onsets),1);
for wave_no = 1 : min([length(onsets), length(peaks)])
    desired_ht = mean(sig([onsets(wave_no), peaks(wave_no)]));
    [~, temp] = min(abs(sig(onsets(wave_no):peaks(wave_no)) - desired_ht));
    mid_amps(wave_no,1) = temp + onsets(wave_no) - 1;
end

end