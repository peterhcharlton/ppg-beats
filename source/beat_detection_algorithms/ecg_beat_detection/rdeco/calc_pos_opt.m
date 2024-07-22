function [pos_opt] = calc_pos_opt(signal,step,maxmin)
% Calculation of optima
%
% Input:    signal      signal
%           step        stepsize according to application
%           maxmin      for finding maxima: maxmin = 1
%                       for finding minima: maxmin = -1
% Output:   pos_opt     positions of the optima
%
% Author(s):    Jonathan Moeyersons       (Jonathan.Moeyersons@esat.kuleuven.be)
%               Sabine Van Huffel         (Sabine.Vanhuffel@esat.kuleuven.be)
%               Carolina Varon            (Carolina.Varon@esat.kuleuven.be)
%
% Version History:
% - 06/05/2019   JM      Initial version
%
% Copyright (c) 2019,  Jonathan Moeyersons, KULeuven-ESAT-STADIUS 
%
% This software is made available for non commercial research purposes only 
% under the GNU General Public License. However, notwithstanding any 
% provision of the GNU General Public License, this software may not be 
% used for commercial purposes without explicit written permission after 
% contacting jonathan.moeyersons@esat.kuleuven.be
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

% Get the length of the signal
sze = length(signal);

%% Select the upward slopes

% Pre-allocate
difference = zeros(sze-1,1);

for ii = 1:sze-step
    % If the amplitude of the selected sample is smaller than the amplitude 
    % of the selected sample plus the step size
    if (signal(ii+step)-signal(ii))*maxmin > 0

        difference(ii,1) = 1;
    else
        difference(ii,1) = 0;
    end
end

% time = 1:length(signal);
% figure; plot(signal); hold on; scatter(time(difference==1),signal(difference==1))

%% Select the optima

% Pre-allocate
pos_opt = zeros(1,sze-step);
no_opt = 0;

for ii = 2:sze-step %step/2+1:sze-step
    real_opt = 1;
    
    % If you are at the end of the upwards slope
    if difference(ii-1,1) == 1 && difference(ii,1) == 0  
        count = 1;
        
        % Check if the upward slope lasts long enough
        while real_opt == 1 && count < step && count < ii
            if difference(ii-count,1) == 1
                real_opt = 1;
            else
                real_opt = 0;
            end
            count = count + 1;
        end
        
        % If the upward slope is long enough
        if real_opt == 1 && count >= 0.75*step
            
            % Select the interval containing the peak
            interval = ii : ii + step; 
            
            % Select the extremum in that interval
            if maxmin == 1
                [~,index] = max(signal(interval));
            else
                [~,index] = min(signal(interval));
            end
            
            % Select the position of the peak
            no_opt = no_opt + 1;
            pos_opt(no_opt) = interval(index);
        end
    end
end

pos_opt = pos_opt(1:no_opt);

% figure
% plot(signal);
% hold on
% plot(pos_opt,signal(pos_opt),'ro')
% title('Optima')