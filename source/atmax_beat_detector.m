function [peaks,onsets] = atmax_beat_detector(sig,fs)
% ATMAX_BEAT_DETECTOR  ATmax PPG beat detector.
%   ATMAX_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
%   using the 'Adaptive Threshold Method (Vmax)' beat detector
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
%   # Reference (algorithm)
%   H. S.﻿Shin et al., 'Adaptive threshold method for the peak detection of photoplethysmographic waveform,' Comput Biol Med Elsevier, vol.39, no.12, pp.1145-52, 2009. <﻿https://doi.org/10.1016/j.compbiomed.2009.10.006>
%   
%   # Reference (implementation)
%   D. Han et al., 'A Real-Time PPG Peak Detection Method for Accurate Determination of Heart Rate during Sinus Rhythm and Cardiac Arrhythmia,' Biosensors, vol.12, no.2, p.82, 2022. <﻿https://doi.org/10.3390/bios12020082>
%   
%   # Author
%   * Dong Han - wrote the 'my_peak_compare_Shin_2009' function.
%   * Peter H. Charlton - did very little, just wrote this wrapper (which detects onsets given the peaks provided by the code)
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Version
%   1.0
%   
%   # Source
%   my_peak_compare_Shin_2009.m, from the PPG_Peak_Detection GitHub Repository (accessed on 3 Mar 2022) at: <https://github.com/Cassey2016/PPG_Peak_Detection>
%   
%   # Licence
%      MIT Licence (see the licence at the top of the 'my_peak_compare_shin_2009' function in the code below).

V_max_flag = 1; % specifies to detect pulse peaks
output = my_peak_compare_Shin_2009(sig,fs,V_max_flag);
peaks = output.PPG_peak_loc_Shin_2009;

peaks = tidy_beats(peaks);

onsets = pulse_onsets_from_peaks(sig, peaks);

end

function [output_Shin_2009] = my_peak_compare_Shin_2009(raw_PPG,fs_PPG,V_max_flag)

% MIT License
% 
% Copyright (c) 2022 Cassey (Dong Han)
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% =========================================================================
% This function is the implementation of this paper:
% Shin, Hang Sik, Chungkeun Lee, and Myoungho Lee. 
% "Adaptive threshold method for the peak detection of 
% photoplethysmographic waveform." 
% Computers in biology and medicine 
% 39.12 (2009): 1145-1152.
% 
% Implemented by: Dong Han, on 02/10/2020.
%
% Please cite our paper if you used this code:
% Han, Dong, Syed K. Bashar, Jesús Lázaro, Fahimeh Mohagheghian, 
% Andrew Peitzsch, Nishat Nishita, Eric Ding, Emily L. Dickson, 
% Danielle DiMezza, Jessica Scott, Cody Whitcomb, Timothy P. Fitzgibbons, 
% David D. McManus, and Ki H. Chon. 2022. 
% "A Real-Time PPG Peak Detection Method for Accurate Determination of 
% Heart Rate during Sinus Rhythm and Cardiac Arrhythmia" 
% Biosensors 12, no. 2: 82. https://doi.org/10.3390/bios12020082 
%
% Please cite our paper if you used our code. Thank you.
% =========================================================================
    debugging_plot_flag = false; % debugging plot. Can be false if don't want to plot anything.
%% Section 2.4 PPG frequency analysis and filtering.

    % (1): high pass >= 0.5 Hz.
    [b, a] = butter(6,[0.5 20]/(fs_PPG/2)); % bandpass filter 0.5-10Hz, changed from 0.5-20 to 0.5-9 Hz at 11/21/2018
    raw_PPG = filtfilt(b, a, raw_PPG); % -> AC component
    raw_PPG = raw_PPG ./ std(raw_PPG); % normalizing data is very important for my peak detection.
    raw_PPG = raw_PPG - mean(raw_PPG);
%% Section 2.5 & 2.6 Peak detection algorithm & Adaptive threshold detection

    % (1): bandpass filtering, no moving average filter or wavelet
    % decomposition.
    filtered_PPG = raw_PPG;
    Fs = fs_PPG;

    % % ===== interpolation to 1kHz of PPG: =====
    % x = 1:length(filtered_PPG);
    % v = filtered_PPG;
    % 
    % upsample_Fs = 250;
    % xq = 1:Fs/upsample_Fs:length(filtered_PPG);
    % vq1 = interp1(x,v,xq);
    % 
    % filtered_PPG = vq1;
    % Fs = upsample_Fs; % upsampled to 1000 Hz.

    % figure
    % plot(x,v,'o',xq,vq1,':.');
    % xlim([0 max(xq)]);
    % title('(Default) Linear Interpolation');

    % (2): V_max
    % slope_k: k-th slope amplitude;
    % s_r: slope changing rate (empirically: V_max = -0.6);
    % V_n_1: previous peak amplitude;
    % std_PPG: standard deviation of entire PPG signal;
    % Fs: sampling frequency.

    filtered_PPG = filtered_PPG(:);
    slope_k = NaN(size(filtered_PPG)); % should be a column vector.
    peak_loc = NaN(size(filtered_PPG)); % the array to store PPG peak index.
    pk_idx = 1; % the counter of peaks.
    %% Section 2.7: Peak Correction
    refractory_period = 0.6 * Fs; % sec * sampling frequency, initial refractory period is 0.6 sec.

    temp_win_left = round(0.15 * Fs); % sec * sampling frequency. This is the search region for local minima or maxima detection. chose 0.15 sec because 0.3 sec == 200 BPM.
    temp_win_right = round(0.15 * Fs);


    if V_max_flag % doing upper peak detection.
        s_r = -0.6; 
    else
        s_r = 0.6;%0.6; % not positive because my signal is zero mean.
        % I need to make all bottom signal positive, so I am moving them up.
    %     move_filter_amp = min(filtered_PPG) * (-1);
    %     filtered_PPG = filtered_PPG + move_filter_amp + std(raw_PPG); % move the lowest value more than zero.
    end


    slope_meet_PPG_flag = false; % mark if the slope meet PPG.
    slope_lower_PPG_flag = false; % mark if slope is lower than PPG, once PPG amp is lower than slope, mark it back.
    prev_slope = NaN; % First, I want to test not decreasing with PPG amplitude.
    if debugging_plot_flag % debugging plot
       figure;
       plot(filtered_PPG);
       hold on;
    end
    for kk = 1:length(filtered_PPG)
        % this is for debugging:
        if kk == 2
            my_stop = 1;
        end
        if kk == 1 % initial the slope value
            if V_max_flag
                slope_k(1,1) = 0.2 * max(filtered_PPG);
                std_PPG = std(filtered_PPG);
            else
                slope_k(1,1) = 0.2 * min(filtered_PPG); % since my signal is zero mean, I start from the negative amp. % I added what I moved.
                std_PPG = -std(filtered_PPG);
            end
    %         std_PPG = std(filtered_PPG);
            V_n_1 = slope_k(1,1);
        else
            if slope_meet_PPG_flag % slope has met PPG before.
                slope_k(kk,1) = filtered_PPG(kk,1);
                if V_max_flag % upper peak detection.
                    if kk < 2 % in the second point of signal
                        turn_point_flag = (slope_k(kk,1) < slope_k(kk-1,1)); % we met local maximum.
                    else
                        turn_point_flag = (slope_k(kk,1) < slope_k(kk-1,1)) & (slope_k(kk - 1,1) > slope_k(kk-2,1)); % we met local maximum.
                    end
                else
                    if kk < 2 % in the second point of signal
                        turn_point_flag = (slope_k(kk,1) > slope_k(kk-1,1)); % we met local minimum.
                    else
                        turn_point_flag = (slope_k(kk,1) > slope_k(kk-1,1)) & (slope_k(kk - 1,1) < slope_k(kk-2,1)); % we met local minimum.
                    end
                end

                if turn_point_flag % there is a turning point.
                    if pk_idx > 1 % not the first peak
                        % check local maxima or minima:
                        if (kk - temp_win_left) < 1
                            temp_left = 1;
                        else
                            temp_left = kk - temp_win_left;
                        end

                        if (kk + temp_win_right) > length(filtered_PPG)
                            temp_right = length(filtered_PPG);
                        else
                            temp_right = kk + temp_win_right;
                        end
                        temp_win = temp_left:temp_right;
                        local_m_check = filtered_PPG(temp_win);
                        if V_max_flag
                            temp_m_idx = find(local_m_check > slope_k(kk - 1,1)); % check if there is another maximum than detected, remember use k-1.
                        else
                            temp_m_idx = find(local_m_check < slope_k(kk - 1,1)); % check if there is another minimum than detected
                        end

                        if isempty(temp_m_idx) % there is no more max or min than this peak
                            if (kk - peak_loc(pk_idx-1,1) > refractory_period) % it is not the first peak, and the second peak is outside refractory period. It should be kk, because I have not assign the peak to the array.
                                peak_loc(pk_idx,1) = kk-1;
                                V_n_1 = filtered_PPG(peak_loc(pk_idx-1,1),1);% previous peak amplitude %slope_k(kk-1,1);
                                % update refractory period:
                                refractory_period = 0.6 * (kk - peak_loc(pk_idx-1,1)); % current index minus peak location. update the refractory peroid before updating the peak counting.
                                pk_idx = pk_idx + 1;

                                % reset slope meet flag:
                                slope_meet_PPG_flag = false;
                                slope_k(kk,1) = slope_k(kk - 1,1) + s_r * ((V_n_1 + std_PPG) / Fs); 

                                % ---- for checking lower slope -------
                                temp_slope_check = s_r * ((V_n_1 + std_PPG) / Fs);
                                if V_max_flag
                                    if temp_slope_check > 0 % upper peaks should be decreasing with negative slope.
                                        temp_slope_check = -s_r * ((V_n_1 + std_PPG) / Fs);%-temp_slope_check;
                                        slope_k(kk,1) = slope_k(kk - 1,1) + temp_slope_check;
                                    end
                                else
                                    if temp_slope_check < 0 % upper peaks should be decreasing with negative slope.
                                        temp_slope_check = -s_r * ((V_n_1 + std_PPG) / Fs);%-temp_slope_check;
                                        slope_k(kk,1) = slope_k(kk - 1,1) + temp_slope_check;
                                    end
                                end
                                % -------------------------------------------
                                if V_max_flag
                                    temp_slope_below_PPG_flag = slope_k(kk,1) < filtered_PPG(kk,1); % upper peak detection, so slope below signal.
                                else
                                    temp_slope_below_PPG_flag = slope_k(kk,1) > filtered_PPG(kk,1); % lower peak detection, so slope above signal.
                                end
                                if temp_slope_below_PPG_flag % if slope is below PPG signal, we will reset slope value to PPG amplitude.
                                    slope_lower_PPG_flag = true; % slope is lower than PPG signal.
                                    prev_slope = slope_k(kk,1); % store the slope value now.
                                    slope_k(kk,1) = filtered_PPG(kk,1);
                                end
                                if debugging_plot_flag % debugging plot
                                    plot(kk,slope_k(kk,1),'r.');
                                end

                            else
                                if (kk - peak_loc(pk_idx-1,1) <= refractory_period) % it is because of the refractory period that cause the no peak. It should be kk, because I have not assign the peak to the array.                      
                                    slope_k(kk,1) = filtered_PPG(kk,1);% from the fig.3(c) in the paper, I see they are using the signal amplitude, not slope. 
                                    % no need to reset slope meet flag, waiting for
                                    % next turning point.
                                    if debugging_plot_flag % debugging plot
                                        plot(kk,slope_k(kk,1),'r.');
                                    end
                                end
                            end
                        else % there are more peaks higher then current kk peak.
                            if debugging_plot_flag % debugging plot
                                plot(kk,slope_k(kk,1),'r.');
                            end
                        end
                    else % the first peak, no need to check refractory period.
                                            % check local maxima or minima:
                        if (kk - temp_win_left) < 1
                            temp_left = 1;
                        else
                            temp_left = kk - temp_win_left;
                        end

                        if (kk + temp_win_right) > length(filtered_PPG)
                            temp_right = length(filtered_PPG);
                        else
                            temp_right = kk + temp_win_right;
                        end
                        temp_win = temp_left:temp_right;
                        local_m_check = filtered_PPG(temp_win);
                        if V_max_flag
                            temp_m_idx = find(local_m_check > slope_k(kk-1,1)); % check if there is another maximum than detected, always detect previous peak.
                        else
                            temp_m_idx = find(local_m_check < slope_k(kk-1,1)); % check if there is another minimum than detected
                        end

                        if isempty(temp_m_idx)
                            peak_loc(pk_idx,1) = kk-1;
                            if pk_idx > 1
                                V_n_1 = filtered_PPG(peak_loc(pk_idx-1,1),1);
                            else
                                V_n_1 = slope_k(kk-1,1);% previous peak amplitude %slope_k(kk-1,1);
                            end
                            pk_idx = pk_idx + 1;

                            % reset slope meet flag:
                            slope_meet_PPG_flag = false;
                            slope_k(kk,1) = slope_k(kk - 1,1) + s_r * ((V_n_1 + std_PPG) / Fs); 
                            % ---- for checking lower slope -------
                            temp_slope_check = s_r * ((V_n_1 + std_PPG) / Fs);
                            if V_max_flag
                                if temp_slope_check > 0 % upper peaks should be decreasing with negative slope.
                                    temp_slope_check = -s_r * ((V_n_1 + std_PPG) / Fs);%-temp_slope_check;
                                    slope_k(kk,1) = slope_k(kk - 1,1) + temp_slope_check;
                                end
                            else
                                if temp_slope_check < 0 % upper peaks should be decreasing with negative slope.
                                    temp_slope_check = -s_r * ((V_n_1 + std_PPG) / Fs);%-temp_slope_check;
                                    slope_k(kk,1) = slope_k(kk - 1,1) + temp_slope_check;
                                end
                            end
                            % -------------------------------------------
                            if V_max_flag
                                temp_slope_below_PPG_flag = slope_k(kk,1) < filtered_PPG(kk,1); % upper peak detection, so slope below signal.
                            else
                                temp_slope_below_PPG_flag = slope_k(kk,1) > filtered_PPG(kk,1); % lower peak detection, so slope above signal.
                            end

                            if temp_slope_below_PPG_flag % if slope is below PPG signal, we will reset slope value to PPG amplitude.
                                slope_k(kk,1) = filtered_PPG(kk,1);
                            end
                                if debugging_plot_flag % debugging plot
                                    plot(kk,slope_k(kk,1),'r.');
                                end
                        else % there are more peaks higher then current kk peak.
                            if debugging_plot_flag % debugging plot
                                plot(kk,slope_k(kk,1),'r.');
                            end
                        end        
                        % no need to calculate refractory period, because there is only one peak, at least two peaks can give this correctly:
                    end
                else
                    % turning point did not meet, so keep decreasing or
                    % increasing the slope.
    %                 slope_k(kk,1) = slope_k(kk - 1,1) + s_r * ((V_n_1 + std_PPG) / Fs); 
                    if debugging_plot_flag % debugging plot
                        plot(kk,slope_k(kk,1),'r.');
                    end
                end
            else % slope has not met PPG before. Keep decresing or increasing according to 'V_max_flag'.
    %             if slope_lower_PPG_flag % if there is a slope lower than PPG before:
    %                 slope_k(kk,1) = prev_slope;
    %             else
                    slope_k(kk,1) = slope_k(kk - 1,1) + s_r * ((V_n_1 + std_PPG) / Fs);
                            % ---- for checking lower slope -------
                            temp_slope_check = s_r * ((V_n_1 + std_PPG) / Fs);
                            if V_max_flag
                                if temp_slope_check > 0 % upper peaks should be decreasing with negative slope.
                                    temp_slope_check = -s_r * ((V_n_1 + std_PPG) / Fs);%-temp_slope_check;
                                    slope_k(kk,1) = slope_k(kk - 1,1) + temp_slope_check;
                                end
                            else
                                if temp_slope_check < 0 % upper peaks should be decreasing with negative slope.
                                    temp_slope_check = -s_r * ((V_n_1 + std_PPG) / Fs);%-temp_slope_check;
                                    slope_k(kk,1) = slope_k(kk - 1,1) + temp_slope_check;
                                end
                            end
                            % -------------------------------------------

    %             end
    %             if slope_k(kk,1) < filtered_PPG(kk,1) % if slope is below PPG signal, we will reset slope value to PPG amplitude.
    %                 slope_lower_PPG_flag = true; % slope is lower than PPG signal.
    %                 prev_slope = slope_k(kk,1); % store the slope value now.
    %                 slope_k(kk,1) = filtered_PPG(kk,1);
    %             elseif slope_k(kk,1) > filtered_PPG(kk,1) % slope is higher.
    %                 slope_lower_PPG_flag = false;
    %                 prev_slope = NaN; % reset the prev value.

    %             end

    %             if slope_lower_PPG_flag ~= 1 % if slope was not lower than PPG.
    %                 % -------------- Check if two lines will meet -----------------
    %                 PPG_x1 = kk - 1;
    %                 PPG_x2 = kk;
    %                 PPG_y1 = filtered_PPG(kk-1,1);
    %                 PPG_y2 = filtered_PPG(kk,1);
    %                 slope = s_r;
    %                 slope_y2 = slope_k(kk,1);
    %                 slope_y1 = slope_k(kk-1,1);
    %                 [meet_x] = my_slope_meet_PPG(PPG_x1,PPG_x2,PPG_y1,PPG_y2,slope,slope_y2,slope_y1);
    % 
    %                 slope_meet_PPG_flag = (ceil(meet_x) == kk);%(slope_k(kk,1) - filtered_PPG(kk,1)) < 0.1; % 0.3 is a testing value. %slope_k(kk,1) == filtered_PPG(kk,1) % slope meets the PPG signal.
    %             end
                if V_max_flag
                    slope_meet_PPG_flag = ((slope_k(kk,1) < filtered_PPG(kk,1)) & slope_k(kk - 1,1) > filtered_PPG(kk - 1,1));
                else
                    slope_meet_PPG_flag = ((slope_k(kk,1) > filtered_PPG(kk,1)) & slope_k(kk - 1,1) < filtered_PPG(kk - 1,1)); % lower peak use inverse amplitude.
                end
                % -------------------------------------------------------------
                % I found I cannot use equal, because the PPG sampling
                % frequency is not so high.
                if slope_meet_PPG_flag
                    slope_k(kk,1) = filtered_PPG(kk,1); % starts from the next index, slope == PPG amplitude.
                else
                    % don't need to do anything.
                    if slope_lower_PPG_flag ~= 1 % there was no slope lower than PPG before.
                        if V_max_flag
                            slope_lower_PPG_flag = ((slope_k(kk,1) < filtered_PPG(kk,1)) & slope_k(kk - 1,1) == filtered_PPG(kk - 1,1)); % beginning part has same amplitude, but the ending part slope is lower.
                        else
                            slope_lower_PPG_flag = ((slope_k(kk,1) > filtered_PPG(kk,1)) & slope_k(kk - 1,1) == filtered_PPG(kk - 1,1)); % lower peak use inverse amplitude.
                        end
                        if slope_lower_PPG_flag
                            prev_slope = slope_k(kk,1); % store the slope value now.
                            slope_k(kk,1) = filtered_PPG(kk,1); % starts from the next index, slope == PPG amplitude.
                        end
                    else % there was slope lower than PPG before.

                        if V_max_flag
                            temp_PPG_below_slope_flag = filtered_PPG(kk,1) < prev_slope; % upper peak detection, so PPG below slope.
                        else
                            temp_PPG_below_slope_flag = filtered_PPG(kk,1) > prev_slope; % lower peak detection, so PPG above slope.
                        end

                        if temp_PPG_below_slope_flag % PPG is lower than prev slope.
                            slope_k(kk,1) = prev_slope; % stop tracking PPG amp.
                            slope_lower_PPG_flag = false; % reset the lower PPG flag.
                            prev_slope = NaN;
                        else
                            slope_k(kk,1) = filtered_PPG(kk,1); % keep tracking PPG amp.
                        end
                    end
                end
                if debugging_plot_flag % debugging plot
                    plot(kk,slope_k(kk,1),'r.');
                end
            end
        end

    end
    % ================== IMPORTANT: clean up NaN value ========================
    peak_loc(isnan(peak_loc)) = []; % remove empty peak loc.
    if V_max_flag % doing upper peak detection.

    else
        % moving signal back.
    %     filtered_PPG = filtered_PPG - move_filter_amp - std(raw_PPG); % move the lowest value more than zero.
    %     slope_k = slope_k - move_filter_amp - std(raw_PPG); % move the slope as well.
    end

    if debugging_plot_flag % debugging plot
        plot(peak_loc,filtered_PPG(peak_loc),'ko');
    end

    if isempty(peak_loc)
        HR_Shin_2009 = 0; % there is no peak location.
        peak_loc = 1;
    else
        HR_Shin_2009 = 60 * Fs ./ diff(peak_loc); % calculate the HR.
    end

    output_Shin_2009 = struct('PPG_peak_loc_Shin_2009',peak_loc,...
        'slope_Shin_2009',slope_k,...
        'filtered_PPG_Shin_2009',filtered_PPG,...
        'HR_Shin_2009',HR_Shin_2009); 
end