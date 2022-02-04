function [peaks, onsets, mid_amps] = detect_ppg_beats(s, options)
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
%   * options - A structure containing the settings used for the processing:
%    - options.beat_detector  - The beat detector to be used
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
switch options.beat_detector
    case 'AMPD'
        peaks = ampd_peak_detector(s.v,s.fs,up);
    case 'IMS'
        [peaks, onsets, ~] = adaptPulseSegment(s.v,s.fs);
        for wave_no = 1 : length(peaks)-1
            [~, temp] = max(s.v(onsets(wave_no):onsets(wave_no+1)));
            peaks(wave_no) = temp+onsets(wave_no)-1;
        end
    case 'ABD'
        peaks = abd_algorithm(s);
    case 'MSPTD'
        [peaks, onsets] = bishop_peak_detector(s);
    case 'Pulses'
        [peaks, ~, peaks_filt, ~] = PPG_pulses_detector(s.v, s.fs);
        peaks = peaks(:);
        % note that I haven't used peaks_filt
    case 'qppg'
        [onsets] = qppg(s.v,s.fs);
        onsets = onsets(:);
    case 'HeartPy'
        peaks = heartpy_pc(s);
    case 'COppg'
        [peaks, onsets] = COppg_peak_detector(s.v, s.fs);
    case 'SPAR3'
        % uses external file
        N = 3;
        [beat_indices_all,cycle_times_all] = SPARcycleTimesPPG(s.v,s.fs,N);
        peaks = nan(length(beat_indices_all)-1,1);
        for beat_no = 1 : length(beat_indices_all)-1
            curr_ind = beat_indices_all(beat_no);
            [~, temp] = max(s.v(curr_ind:beat_indices_all(beat_no+1)));
            peaks(beat_no) = curr_ind+temp-1;
        end
    case 'SPAR7'
        % uses external file
        N = 7;
        [beat_indices_all,cycle_times_all] = SPARcycleTimesPPG(s.v,s.fs,N);        
        peaks = nan(length(beat_indices_all)-1,1);
        for beat_no = 1 : length(beat_indices_all)-1
            curr_ind = beat_indices_all(beat_no);
            [~, temp] = max(s.v(curr_ind:beat_indices_all(beat_no+1)));
            peaks(beat_no) = curr_ind+temp-1;
        end
    case 'ERMA'
        peaks = erma(s);
    case 'PWD'
        [onsets, peaks] = pwd(s);
    case 'PDA'
        peaks = pda(s);
    case 'WFD'
        peaks = wfd(s);
end

% tidy up peaks and onsets
if exist('onsets','var') && exist('peaks','var')
    
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
end

% Add onsets if required
if ~exist('onsets','var')
    onsets = nan(length(peaks)-1,1);
    for wave_no = 1 : length(peaks)-1
        [~, temp] = min(s.v(peaks(wave_no):peaks(wave_no+1)));
        onsets(wave_no) = temp + peaks(wave_no) - 1;
    end
    peaks = peaks(2:end); % so that there are the same number of peaks and onsets
end

% Add peaks if required
if ~exist('peaks','var')
    peaks = nan(length(onsets)-1,1);
    for wave_no = 1 : length(onsets)-1
        [~, temp] = max(s.v(onsets(wave_no):onsets(wave_no+1)));
        peaks(wave_no,1) = temp + onsets(wave_no) - 1;
    end
    onsets = onsets(1:end-1); % so that there are the same number of peaks and onsets
end

% Add mid-amplitude point if required
if ~exist('mid_amps','var')
    mid_amps = nan(length(onsets),1);
    for wave_no = 1 : length(onsets)
        desired_ht = mean(s.v([onsets(wave_no), peaks(wave_no)]));
        [~, temp] = min(abs(s.v(onsets(wave_no):peaks(wave_no)) - desired_ht));
        mid_amps(wave_no,1) = temp + onsets(wave_no) - 1;
    end
end

end
