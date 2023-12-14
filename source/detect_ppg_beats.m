function [peaks, onsets, mid_amps, t_taken] = detect_ppg_beats(s, beat_detector, do_timing)
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
%   * beat_detector  - a string specifying the beat detector algorithm to be used
%
%   * do_timing - a logical indicating whether or not to time how long it takes to run the beat detector algorithm
%   
%   # Outputs
%   * peaks : indices of pulse peaks
%   * onsets : indices of pulse onsets
%   * mid_amps : indices of mid-points on the upslope between onsets and peaks (defined as the mid-amplitude)
%   * t_taken : time taken (in secs) to run the beat detector algorithm
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Author
%   Peter H. Charlton, University of Cambridge, January 2023.
%   
%   # Version
%   1.2
%   
%   # License - GPL-3.0
%      Copyright (c) 2022 Peter H. Charlton
%      This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%      This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%      You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

%% Setup
if nargin<3
    do_timing = 0;
end

%% Detect peaks and onsets in pulse signal
if do_timing, tic, end
eval(['[peaks, onsets] = ' lower(beat_detector) '_beat_detector(s.v,s.fs);']);
if do_timing, t_taken = toc; else, t_taken = nan; end

%% Tidy up peaks and onsets
[peaks, onsets] = tidy_peaks_and_onsets(s.v, peaks,onsets);

%% Add mid-amplitude points
mid_amps = calc_mid_amp_points(s.v, peaks, onsets);

end

function [peaks, onsets] = tidy_peaks_and_onsets(sig, peaks, onsets)

orig_peaks = peaks; orig_onsets = onsets;

% Tidy up peaks and onsets to make them conform to the following rules:
% (i) No two points at the same time
% (ii) At least one local minimum between consecutive peaks
% (iii) At least one local maximum between consecutive onsets
% (iv) Alternates between onsets and peaks
% (v) Starts with onset, and ends with peak
% (vi) Same number of peaks and onsets

%% (i) No two points at the same time
[peaks,onsets] = remove_repeated_peaks_and_onsets(peaks, onsets);

%% (ii) At least one local minimum between consecutive peaks
peaks = ensure_at_least_one_extremum_between_other_type_of_extremum(sig,peaks, 'pk');

%% (iii) At least one local maximum between consecutive onsets
onsets = ensure_at_least_one_extremum_between_other_type_of_extremum(sig,onsets, 'on');

%% (iv) Alternates between onsets and peaks
% If there are two consecutive peaks, then insert an onset between them
[onsets, peaks] = insert_extremum_between_other_type_of_extremum(sig, onsets, peaks, 'pk');
% If there are two consecutive onsets, then insert a peak between them
[peaks, onsets] = insert_extremum_between_other_type_of_extremum(sig, peaks, onsets, 'on');

%% (v) Starts with onset, and ends with peak
[peaks, onsets] = ensure_starts_with_onset_ends_with_peak(peaks,onsets);

%% (vi) same number of peaks and onsets
[peaks, onsets] = ensure_same_no_peaks_onsets(peaks,onsets);

%% check that the peaks and onsets now satisfy the requirements
failed_test = do_test(peaks,onsets);

end

function failed_test = do_test(peaks,onsets, test_types)

fprintf('\n --- Doing tests ---')
failed_test = false;
% Check whether peaks and onsets meet the following rules:
if nargin<3
    test_types = {'duplicates', 'alternates', 'start_end', 'same_no'};
end

% (i) No two points at the same time (duplicates)
if sum(contains(test_types, 'duplicates'))
    all_extrema = [peaks;onsets];
    [~, I] = unique(all_extrema, 'first'); % from: https://uk.mathworks.com/matlabcentral/answers/13149-finding-duplicates#answer_17970
    x = 1:length(all_extrema);
    x(I) = [];
    if ~isempty(x)
        fprintf('\n - Failed on test: there are duplicate points')
        failed_test = true;
    end
end

% (iv) Alternates between onsets and peaks (alternates)
if sum(contains(test_types, 'alternates'))
    peak_log = [true(size(peaks)); false(size(onsets))];
    all_extrema = [peaks;onsets];
    [~, order] = sort(all_extrema);
    peak_log_ordered = peak_log(order);
    if sum(abs(diff(peak_log_ordered))~=1)
        fprintf('\n - Failed on test: does not alternate between onsets and peaks')
        failed_test = true;
    end
end

% (v) Starts with onset, and ends with peak (start_end)
if sum(contains(test_types, 'start_end'))
    if min(onsets) > min(peaks)
        fprintf('\n - Failed on test: does not start with onset')
        failed_test = true;
    end
    if max(onsets) > max(peaks)
        fprintf('\n - Failed on test: does not end with peak')
        failed_test = true;
    end
end

% (vi) Same number of peaks and onsets (same_no)
if sum(contains(test_types, 'same_no'))
    if length(onsets) ~= length(peaks)
        fprintf('\n - Failed on test: does not have same number of onsets and peaks')
        failed_test = true;
    end
end

if ~failed_test
    fprintf(' Passed tests ---')
else
    fprintf(' Failed tests ---')
end

end

function [peaks, onsets] = ensure_same_no_peaks_onsets(peaks,onsets)

% NB: This doesn't quite ensure the same no of peaks and onsets, it only does it for a specific condition

% if no peaks (or onsets) were detected, then don't output any indices for either
if isempty(peaks)
    onsets = [];
end
if isempty(onsets)
    peaks = [];
end

end

function [peaks, onsets] = ensure_starts_with_onset_ends_with_peak(peaks,onsets)
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

end

function [extrema, other_extrema] = insert_extremum_between_other_type_of_extremum(sig, extrema, other_extrema, other_extrema_type)

% If there are two consecutive extrema of one type (known as other extrema), then insert an extremum of the second type (known as extrema) between them
other_extrema_log = [true(size(other_extrema)); false(size(extrema))];
[els,order] = sort([other_extrema;extrema]);
other_extrema_log = other_extrema_log(order);
bad_els = find(diff(other_extrema_log)==0 & other_extrema_log(1:end-1)); % repeated other_extrema
if ~isempty(bad_els)  % if there is a repeated other extrema
    for bad_el_no = 1 : length(bad_els)   % cycle through each repeated other extrema
        curr_other_extrema = [els(bad_els(bad_el_no)),els(bad_els(bad_el_no)+1)];
        bw_to_remove = linspace(sig(curr_other_extrema(1)),sig(curr_other_extrema(2)), curr_other_extrema(2)-curr_other_extrema(1)+1 );
        if strcmp(other_extrema_type, 'pk')
            [~,temp]=min(sig(curr_other_extrema(1):curr_other_extrema(2)) - bw_to_remove(:));
        else
            [~,temp]=max(sig(curr_other_extrema(1):curr_other_extrema(2)) - bw_to_remove(:));
        end
        % check this hasn't just detected one of the other_extrema (which can happen with strong baseline wander)
        if temp==1 || temp==(curr_other_extrema(2)-curr_other_extrema(1)+1)
            % then just remove the first peak
            other_extrema(other_extrema==curr_other_extrema(1)) = [];
        else
            curr_new_extrema = curr_other_extrema(1) -1 + temp;
            extrema = sort([extrema; curr_new_extrema]);
        end
        
    end
end

end

function other_extrema = ensure_at_least_one_extremum_between_other_type_of_extremum(sig, other_extrema, other_extrema_type)

% If there are two peaks (or onsets) without a local minimum (or maximum) between them, then remove the lower (or higher) one

if strcmp(other_extrema_type, 'pk')
    extrema = islocalmin(sig);
else
    extrema = islocalmax(sig);
end
finished = false;
while ~finished
    els_to_remove = [];
    for other_extrema_no = 1 : length(other_extrema)-1
        rel_els = (other_extrema(other_extrema_no)+1) : (other_extrema(other_extrema_no+1)-1);
        if sum(extrema(rel_els)) == 0
            other_extrema_vals = sig(other_extrema([other_extrema_no,other_extrema_no+1]));
            if strcmp(other_extrema_type, 'pk')
                [~, el_to_remove] = min(other_extrema_vals);
            else
                [~, el_to_remove] = max(other_extrema_vals);
            end
            els_to_remove(end+1) = other_extrema_no + el_to_remove - 1;
        end
    end
    other_extrema(els_to_remove) = [];
    if isempty(els_to_remove)
        finished = true;
    end
end

% % If there are two onsets without a local maximum between them, then remove the higher one
% local_max = islocalmax(sig);
% finished = false;
% while ~finished
%     els_to_remove = [];
%     for onset_no = 1 : length(onsets)-1
%         rel_els = (onsets(onset_no)+1) : (onsets(onset_no+1)-1);
%         if sum(local_max(rel_els)) == 0
%             onset_vals = sig(onsets([onset_no,onset_no+1]));
%             [~, el_to_remove] = max(onset_vals);
%             els_to_remove(end+1) = onset_no + el_to_remove - 1;
%         end
%     end 
%     onsets(els_to_remove) = [];
%     if isempty(els_to_remove)
%         finished = true;
%     end
% end
% clear local_max finished els_to_remove onset_no rel_els onset_vals el_to_remove


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

function [peaks,onsets] = remove_repeated_peaks_and_onsets(peaks, onsets)

% remove any repeated peaks (or onsets)
peaks = unique(peaks);
onsets = unique(onsets);

% If there is a peak and onset at the same index, then remove them both
repeated_vals = intersect(peaks,onsets);
peaks = setxor(peaks, repeated_vals);
onsets = setxor(onsets, repeated_vals);
clear repeated_vals

end