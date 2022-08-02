function [peaks, onsets] = erma_beat_detector(sig,fs)
% ERMA_BEAT_DETECTOR  ERMA PPG beat detector.
%   ERMA_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
%   using the 'Event-Related Moving Averages' beat detector
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
%   M. Elgendi et al., 'Systolic peak detection in acceleration photoplethysmograms measured from emergency responders in tropical conditions,' PLoS ONE, vol. 8, no. 10, pp. 1-11, 2013. <https://doi.org/10.1371/journal.pone.0076585>
%   
%   # Author
%   * Elisa Mejía Mejía - wrote the code: 'erma'
%   * Peter H. Charlton - did very little, just wrote this wrapper and made slight edits (see 'erma')
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Version
%   1.0
%   
%   # MIT License
%      Copyright (c) 2022 Elisa Mejía Mejía and Peter H. Charlton
%      Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
%      The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
%      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

s.v = sig;
s.fs = fs;
peaks = erma(s);

peaks = tidy_beats(peaks);

onsets = pulse_onsets_from_peaks(sig, peaks);

end

function peaks = erma(s)

%% Detects inter-beat intervals using D2Max
% Citation: Elgendi M, Norton I, Brearley M, Abbott D, Schuurmans D (2013) Systolic Peak Detection in Acceleration
% Photoplethysmograms Measured from Emergency Responders in Tropical Conditions. PLoS ONE 8(10): e76585.
% doi:10.1371/journal.pone.0076585
%
% Inputs:   PPG signal, ppg
%           sampling rate, fs
%           visualization option, v (1: plot, 0: don't plot)
% Outputs:  position of the starting points of inter-beat intervals, in number of samples, ibis
%
% Developed by: Elisa Mejía Mejía
%               City, University of London
%
% Adapted by: Peter H. Charlton
%
% Version:      1.1 -   Nov, 2021

%% setup 
% (inserted by PHC)
ppg = s.v;
fs = s.fs;
clear s
v = 0; % set to 1 to enable visualiation

%% Bandpass filter
% if length(ppg) < 4098                                   % Removed by PHC
%     ppg = [ppg; zeros(4098 - length(ppg) + 1,1)];
% end
[b,a] = butter(2,[0.5, 8]/(fs/2),'bandpass');   % Designs the bandpass, butterworth filter
s = filtfilt(b,a,ppg);                          % Applies the zero-phase filter
z = s;                                          % Creates a new signal based on the filtered signal
z(z < 0) = 0;                                   % Clips the signal

if v == 1
    color = [   0, 0.4470, 0.7410; ...
        0.8500, 0.3250, 0.0980; ...
        0.9290, 0.6940, 0.1250; ...
        0.4940, 0.1840, 0.5560; ...
        0.4660, 0.6740, 0.1880; ...
        0.3010, 0.7450, 0.9330; ...
        0.6350, 0.0780, 0.1840; ...
        0.75, 0, 0.75];
    
    figure('position',[240 140 1200 625],'name','D2Max');
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
    plot(s,'k','linewidth',1.2);
    hold on;
    plot(z,'color',color(1,:),'linewidth',1.2);
    xlim([1 length(s)]);
    title('D2Max');
    xlabel('Samples');
    ylabel('Filtered PPG');
    grid on;
    box off;
    linkaxes(ax,'x');
end

%% Squaring
y = z.^2;                                       % Squares the filtered, clipped signal

if v == 1
    ax(2) = subplot(2,1,2);
    plot(y,'color',color(2,:),'linewidth',1.2);
    linkaxes(ax,'x');
end

%% Generating blocks of interest
w1 = (111e-3)*fs;                               % Determines the order of a first moving average
w1 = 2*floor(w1/2) + 1;                         % Rounds the order to the nearest odd integer
b = (1/w1)*ones(w1,1);                          % Designs the first moving average filter
ma_peak = filtfilt(b,1,y);                      % Applies the first moving average filter

w2 = (667e-3)*fs;                               % Determines the order of a second moving average
w2 = 2*floor(w2/2) + 1;                         % Rounds the order to the nearest odd integer
b = (1/w2)*ones(w2,1);                          % Designs the second moving average filter
ma_beat = filtfilt(b,1,y);                      % Applies the second moving average filter

if v == 1
    ax(2) = subplot(2,1,2);
    plot(ma_peak,'--','color',color(3,:),'linewidth',1.2);
    plot(ma_beat,'--','color',color(4,:),'linewidth',1.2);
    linkaxes(ax,'x');
end

%% Thresholding
alpha = 0.02*mean(y);                           % Determines offset level
th1 = ma_beat + alpha;                          % Determines a first threshold based on crossing of ma_beat
boi = ma_peak > th1;                            % Determines blocks of interest according to threshold 1

if v == 1
    ax(2) = subplot(2,1,2);
    plot(boi,'color',color(5,:),'linewidth',1.2);
    linkaxes(ax,'x');
end

th2 = w1;                                       % Determines a second threshold based on length of the
% blocks of interest
pos_blocks_init = find(diff(boi) > 0);          % Determines the location of the starting points for blocks
pos_blocks_init = pos_blocks_init + 1;          % Adds one for the delay
pos_blocks_end = find(diff(boi) < 0);           % Determines the location of the ending points for blocks
pos_blocks_end = pos_blocks_end + 1;            % Adds one for the delay
if pos_blocks_init(1) > pos_blocks_end(1)       % Verifies if the first block does not have a starting point
    pos_blocks_init = [1, pos_blocks_init];     % Adds a initial point for the first block
end
if pos_blocks_init(end) > pos_blocks_end(end)   % Verifies if the last block does not have an ending point
    pos_blocks_end = [pos_blocks_end, length(y)];
    % Adds an ending point for the last block
end

len_blocks = zeros(size(pos_blocks_init));      % Initializes variable
ibis = zeros(size(pos_blocks_init));            % Initializes variable
for i = 1:length(pos_blocks_init)               % Iterates through the blocks
    ind = find(pos_blocks_end > pos_blocks_init(i),1);
    % Determines the first ending position occurring after
    % the i-th initial position
    len_blocks(i) = pos_blocks_end(ind) - pos_blocks_init(i);
    % Measures the size of the block
    if len_blocks(i) >= th2                     % Verifies if the block is long enough
        [~,max_block] = max(ppg(pos_blocks_init(i):pos_blocks_end(ind)));
        % Detects the maximum in the signal within the block
        ibis(i) = max_block + pos_blocks_init(i) - 1;
        % Adds the offset and stores the results
        
    end
end
ind = find(len_blocks < th2);                   % Finds the blocks that are not long enough
if ~isempty(ind)                                % Checks if at least one block needs to be deleted
    for i = 1:length(ind)                       % Iterates through the blocks to be deleted
        boi(pos_blocks_init(i):pos_blocks_end(i)) = 0;
    end
end
ibis(ibis == 0) = [];                           % Removes the zeros from the array

if v == 1
    ax(1) = subplot(2,1,1);
    plot(ibis,ppg(ibis),'ok','MarkerFaceColor',color(1,:));
    ax(2) = subplot(2,1,2);
    plot(boi,'--','color',color(6,:),'linewidth',1.2);
    plot(ibis,s(ibis),'ok','MarkerFaceColor',color(7,:));
    linkaxes(ax,'x');
end

%% End matter
% (inserted by PHC)
peaks = ibis;

end
