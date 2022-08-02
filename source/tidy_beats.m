function beat_indices = tidy_beats(beat_indices)
% TIDY_BEATS  Cleans up beat detections.
%   TIDY_BEATS makes the vector of beat indices into a column vector of
%   unique values.
%   
%   # Inputs
%   
%   * beat_indices  : a vector containing the indices of beat detections
%   
%   # Outputs
%   * beat_indices  : a cleaned up vector of indices of beat detections
%   
%   # Author
%   Peter H. Charlton, February 2022
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

% - ensure it's a column vector
beat_indices = beat_indices(:);

% - only keep unique values
beat_indices = unique(beat_indices);

end