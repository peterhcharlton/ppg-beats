function [peaks, onsets] = pda_beat_detector(sig,fs)
% PDA_BEAT_DETECTOR  PDA PPG beat detector.
%   PDA_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
%   using the 'Peak Detection Algorithm' beat detector
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
%   E. J. Arguello Prada and R. D. Serna Maldonado, 'A novel and low-complexity peak detection algorithm for heart rate estimation from low-amplitude photoplethysmographic (PPG) signals,' Journal of Medical Engineering and Technology, vol. 42, no. 8, pp. 569-577, 2018. <https://doi.org/10.1080/03091902.2019.1572237>
%   
%   # Author
%   * Elisa Mejía Mejía - wrote the code: 'upslopes'
%   * Peter H. Charlton - did very little, just wrote this wrapper
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Version
%   1.0
%   
%   # License - MIT
%      Copyright (c) 2022 Elisa Mejía Mejía and Peter H. Charlton
%      Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
%      The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
%      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

v = 0; % set to 1 to enable visualiation
peaks = upslopes(sig,v);

peaks = tidy_beats(peaks);

onsets = pulse_onsets_from_peaks(sig, peaks);

end

function ibis = upslopes(ppg,v)
%% Detects inter-beat intervals using the upslopes
% Citation: Argüello Prada EJ, Serna Maldonado RD. (2018 ) A novel and low-complexity peak detection algorithm for 
% heart rate estimation from low-amplitude photoplethysmographic (PPG) signals. J Med Eng Technol 42(8):569-577. 
% doi: 10.1080/03091902.2019.1572237. Epub 2019 Mar 28. PMID: 30920315.
%
% Inputs:   PPG signal, ppg
%           visualization option, v (1: plot, 0: don't plot)
% Outputs:  position of the starting points of inter-beat intervals, in number of samples, ibis
% 
% Developed by: Elisa Mejía Mejía
%               City, University of London
% Version:      1.0 -   May, 2021

    %% Detection of peaks 
    if v == 1
        color = [   0, 0.4470, 0.7410; ...
                    0.8500, 0.3250, 0.0980; ...
                    0.9290, 0.6940, 0.1250; ...
                    0.4940, 0.1840, 0.5560; ...
                    0.4660, 0.6740, 0.1880; ...
                    0.3010, 0.7450, 0.9330; ...
                    0.6350, 0.0780, 0.1840; ...
                    0.75, 0, 0.75];
                
        figure('position',[240 140 1200 625],'name','Upslopes');
        ax(1) = subplot(2,1,1);
        plot(ppg,'k','linewidth',1.2);
        hold on;
        xlim([1 length(ppg)]);
        title('Peak detection');
        xlabel('Samples');
        ylabel('PPG');
        grid on;
        box off;
        
        ax(2) = subplot(2,1,2);
        plot(ppg,'k','linewidth',1.2);
        hold on;
        xlim([1 length(ppg)]);
        title('Upslopes');
        xlabel('Samples');
        ylabel('PPG');
        grid on;
        box off;  
        linkaxes(ax,'x');   
    end
    
    th = 6;                                     % Initializes threshold
    pks = [];                                   % Initializes variable  
    pos_peak = [];                              % Initializes variable
    pos_peak_b = 0;                             % Initializes variable
    n_pos_peak = 0;                             % Initializes counter
    n_up = 0;                                   % Initializes counter
    for i = 2:length(ppg)                       % Iterates through signal
        if ppg(i) > ppg(i - 1)                  % Verifies if it is the upslope
            n_up = n_up + 1;                    % Adds one to the coutner
            
            if v == 1
                ax(2) = subplot(2,1,2);
                plot(i,ppg(i),'.','color',color(5,:));
                linkaxes(ax,'x');
            end
        else
            if n_up >= th                       % Checks if the number of data in the upslope is higher than the
                                                % threshold
                pos_peak = [pos_peak; i];       % Adds the position to the possible peaks 
                pos_peak_b = 1;                 % Sets the value as 1 for the detection of a new possible peak
                n_pos_peak = n_pos_peak + 1;    % Refreshes counter
                n_up_pre = n_up;                % Stores the previous value of n_up
                
                if v == 1
                    ax(2) = subplot(2,1,2);
                    plot(pos_peak(n_pos_peak),ppg(pos_peak(n_pos_peak)),'ok','MarkerFaceColor',color(1,:));
                    linkaxes(ax,'x');
                end
            else
                if pos_peak_b == 1              % Verifies if a peak has been found
                    if ppg(i - 1) > ppg(pos_peak(n_pos_peak))
                                                % Verifies if the previous sample is higher than
                                                % the one selected as peak
                        pos_peak(n_pos_peak) = i - 1;
                                                % Updates the value of the position of the possible
                                                % peak
                        if v == 1
                            ax(2) = subplot(2,1,2);
                            plot(pos_peak(n_pos_peak),ppg(pos_peak(n_pos_peak)),'ok','MarkerFaceColor',color(2,:));
                            linkaxes(ax,'x');
                        end
                    else
                        pks = [pks; pos_peak(n_pos_peak)];
                                                % Refreshes array of peaks
                        if v == 1
                            ax(2) = subplot(2,1,2);
                            plot(pks(end),ppg(pks(end)),'ok','MarkerFaceColor',color(3,:));
                            linkaxes(ax,'x');
                        end
                    end
                    th = 0.6*n_up_pre;          % Updates value of threshold
                    pos_peak_b = 0;             % Updates value of variable
                end
            end
            n_up = 0;                           % Resets counter
        end
    end
    
    if v == 1
        ax(1) = subplot(2,1,1);
        plot(pks,ppg(pks),'ok','MarkerFaceColor',color(1,:));
        linkaxes(ax,'x');
    end
    
    ibis = pks;                                 % Sets the output
    
end