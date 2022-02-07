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
%   * peaks : ...TBC...
%   * onsets : ...TBC...
%   * mid_amps : ...TBC...
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Author
%   Peter H. Charlton, University of Cambridge, February 2022.
%   
%   # Version
%   0.1, and is still in development.
%   
%   # Licence
%      This file is part of PPG-beats.
%      PPG-beats is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%      PPG-beats is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%      You should have received a copy of the GNU General Public License along with PPG-beats. If not, see <https://www.gnu.org/licenses/>.

%% Detect peaks in pulse signal

eval(['[peaks, onsets] = ' lower(beat_detector) '_beat_detector(s.v,s.fs);']);

%% Tidy up peaks and onsets

% If there is a peak and onset at the same index, then remove them
peak_log = [true(size(peaks)); false(size(onsets))];
[els,order] = sort([peaks;onsets]);
peak_log = peak_log(order);
bad_els = find(diff(els)==0 & peak_log(1:end-1)~=peak_log(2:end)); % position of peak and onset at the same index
[peaks_to_remove, onsets_to_remove] = deal([]);
if ~isempty(bad_els)  % if there is a repeated peak
    for bad_el_no = 1 : length(bad_els)   % cycle through each repeated peak
        curr_index = els(bad_els(bad_el_no));
        peaks = peaks(peaks~=curr_index);
        onsets = onsets(onsets~=curr_index);
    end
    
    peaks = peaks(~peaks_to_remove);
end
clear new_troughs temp curr_pks bad_el_no bad_els peak_log order els

% If there are two consecutive peaks, then insert a trough between them
peak_log = [true(size(peaks)); false(size(onsets))];
[els,order] = sort([peaks;onsets]);
peak_log = peak_log(order);
bad_els = find(diff(peak_log)==0 & peak_log(1:end-1)); % repeated peaks
if ~isempty(bad_els)  % if there is a repeated peak
    new_troughs = nan(length(bad_els),1);
    for bad_el_no = 1 : length(bad_els)   % cycle through each repeated peak
        curr_pks = [els(bad_els(bad_el_no)),els(bad_els(bad_el_no)+1)];
        [~,temp]=min(s.v(curr_pks(1):curr_pks(2)));
        new_troughs(bad_el_no,1) = curr_pks(1) -1 + temp;
    end
    onsets = sort([onsets; new_troughs]);
end
clear new_troughs temp curr_pks bad_el_no bad_els peak_log order els

% If there are two consecutive troughs, then insert a peak between them
trough_log = [false(size(peaks)); true(size(onsets))];
[els,order] = sort([peaks;onsets]);
trough_log = trough_log(order);
bad_els = find(diff(trough_log)==0 & trough_log(1:end-1)); % repeated troughs
if ~isempty(bad_els)  % if there is a repeated peak
    new_peaks = nan(length(bad_els),1);
    for bad_el_no = 1 : length(bad_els)   % cycle through each repeated peak
        curr_trs = [els(bad_els(bad_el_no)),els(bad_els(bad_el_no)+1)];
        [~,temp]=max(s.v(curr_trs(1):curr_trs(2)));
        new_peaks(bad_el_no,1) = curr_trs(1) -1 + temp;
    end
    peaks = sort([peaks; new_peaks]);
end
clear new_peaks temp curr_trs bad_el_no bad_els trough_log order els

% Make sure that the first onset is before the first peak, and the last peak is after the last onset
if ~isempty(onsets) && ~isempty(peaks)
    if onsets(1)>peaks(1)
        peaks(1) = [];
    end
    if peaks(end)<onsets(end)
        onsets(end) = [];
    end
end

% Add mid-amplitude point if required
if ~exist('mid_amps','var')
    mid_amps = nan(length(onsets),1);
    for wave_no = 1 : min([length(onsets), length(peaks)])
        desired_ht = mean(s.v([onsets(wave_no), peaks(wave_no)]));
        [~, temp] = min(abs(s.v(onsets(wave_no):peaks(wave_no)) - desired_ht));
        mid_amps(wave_no,1) = temp + onsets(wave_no) - 1;
    end
end

end
