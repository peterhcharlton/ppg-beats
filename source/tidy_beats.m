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
%   0.1, and is still in development.
%   
%   # Licence
%      This is available under the GNU General Public License v2.0: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html

% - ensure it's a column vector
beat_indices = beat_indices(:);

% - only keep unique values
beat_indices = unique(beat_indices);

end