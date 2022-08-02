function peaks = pulse_peaks_from_onsets(sig, onsets)
% PULSE_PEAKS_FROM_ONSETS  Identifies pulse onsets.
%   PULSE_PEAKS_FROM_ONSETS detects pulse peaks in a photoplethysmogram 
%   (PPG) signal from the locations of the pulse onsets
%   
%   # Inputs
%   
%   * sig : a vector of PPG values
%   * onsets : indices of detected pulse onsets
%   
%   # Outputs
%   * peaks : indices of detected pulse peaks
%   
%   # Author
%   Peter H. Charlton, University of Cambridge, February 2022
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Version
%   1.0
%   
%   # License - GPL-3.0
%      Copyright (c) 2022 Peter H. Charlton
%      This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%      This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%      You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

peaks = nan(length(onsets)-1,1);
for wave_no = 1 : length(onsets)-1
    [~, temp] = max(sig(onsets(wave_no):onsets(wave_no+1)));
    peaks(wave_no) = temp + onsets(wave_no) - 1;
end

end