function [peaks, onsets] = wfd_beat_detector(sig,fs)
% WFD_BEAT_DETECTOR  WFD PPG beat detector.
%   WFD_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
%   using the 'Wavelet Foot Delineation' beat detector
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
%   N. J. Conn and D. A. Borkholder, 'Wavelet based photoplethysmogram foot delineation for heart rate variability applications,' in IEEE Signal Processing in Medicine and Biology Symposium. IEEE, 2013. <https://doi.org/10.1109/SPMB.2013.6736782>
%   
%   # Author
%   * Elisa Mejía Mejía - wrote the code: 'wavelet'
%   * Peter H. Charlton - did very little, just wrote this wrapper and made slight edits (see 'erma')
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
peaks = wavelet(sig,fs,v);

peaks = tidy_beats(peaks);

onsets = pulse_onsets_from_peaks(sig, peaks);
peaks = pulse_peaks_from_onsets(sig, onsets);  % this was needed because the original beat detections were on upslopes rather than peaks

end

function ibis = wavelet(ppg,fs,v)
%% Detects inter-beat intervals using a Wavelet approach
% Citation: Conn NJ, Borkholder DA (2013) Wavelet based photoplethysmogram foot delineation for heart rate 
% variability applications. IEEE Signal Processing in Medicine and Biology Symposium (SPMB) 2013:1-5. 
% doi: 10.1109/SPMB.2013.6736782.
%
% Inputs:   PPG signal, ppg
%           sampling rate, fs
%           visualization option, v (1: plot, 0: don't plot)
% Outputs:  position of the starting points of inter-beat intervals, in number of samples, ibis
% 
% Developed by: Elisa Mejía Mejía
%               City, University of London
% Version:      1.0 -   May, 2021

    %% Beat segmentation
    b = fir1(round(fs/10),[0.5 8]/(fs/2),'bandpass');   % Designs the filter
    ppgf = filtfilt(b,1,ppg);                           % Filters the signal   
    
    if v == 1
        color = [   0, 0.4470, 0.7410; ...
                    0.8500, 0.3250, 0.0980; ...
                    0.9290, 0.6940, 0.1250; ...
                    0.4940, 0.1840, 0.5560; ...
                    0.4660, 0.6740, 0.1880; ...
                    0.3010, 0.7450, 0.9330; ...
                    0.6350, 0.0780, 0.1840; ...
                    0.75, 0, 0.75; ...
                    0.5, 0.5, 0.5];
                
        figure('position',[240 140 1200 625],'name','Wavelet');
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
        yyaxis left;
        plot(ppgf,'-k','linewidth',1.2);
        hold on;
        xlim([1 length(ppg)]);
        title('Wavelet');
        xlabel('Samples');
        ylabel('PPG');
        grid on;
        box off;  
        
        linkaxes(ax,'x'); 
    end
    
    fsi = 250;                                          % Defines the interpolation frequency
    ti = (0:length(ppgf) - 1)/fs;                       % Generates the original time array
    tf = ti(1):1/fsi:ti(end);                           % Generates the time array for the spline
    ppgss = spline(ti,ppgf,tf);                         % Subsamples the PPG signal to 250 Hz
    
    [~,q] = wavelet_filters(length(ppgss),fsi,0);       % Generates the Wavelet filters at 250 Hz
    w = wavelet_decomposition(ppgss,q,0);               % Applies the Wavelet decomposition to the 
                                                        % subsampled signal
    ppgw = w{5};                                        % Takes the fifth Wavelet scale
    ppgw = spline(tf,ppgw,ti);                          % Restores sampling rate
    
    if v == 1
        ax(2) = subplot(2,1,2);
        yyaxis left;
        plot(ppgw,'-','color',color(2,:),'linewidth',1.2);
        xlim([1 length(ppgw)]);
        
        linkaxes(ax,'x');  
    end
    
    ppgrs = ppgw;                                       % Creates a copy of the signal
    ppgrs(ppgrs < 0) = 0;                               % Rectifies the signal
    ppgrs = ppgrs.^2;                                   % Squares the rectified signal
    
    if v == 1
        ax(2) = subplot(2,1,2);
        yyaxis left;
        plot(ppgrs,'-','color',color(3,:),'linewidth',1.2);
        xlim([1 length(ppgrs)]); 
        
        linkaxes(ax,'x');  
    end
    
    b = fir1(round(length(ppgrs)/10),1/(fs/2),'low');   % Designs a Hamming FIR filter
    th = filtfilt(b,1,ppgrs);                           % Applies the Hamming FIR filter
    
    if v == 1
        ax(2) = subplot(2,1,2);
        yyaxis left;
        plot(th,'--','color',color(4,:),'linewidth',1.2);
        xlim([1 length(th)]);
        
        linkaxes(ax,'x');  
    end
    
    zc = [0 diff(ppgrs >= th)];                         % Finds the points where the Wavelet signal crosses the
                                                        % threshold
    zc_up = find(zc == 1);                              % Finds the left boundaries of the areas of interest
    zc_down = find(zc == -1);                           % Finds the right boundaries of the areas of interest
    if zc_up(1) > zc_down(1)                            % Verifies if the first left boundary is higher than
                                                        % the first right boundary
        zc_down(1) = [];                                % Removes the first right boundary
    end
    if zc_up(end) > zc_down(end)                        % Verifies if the last left boundary is higher than
                                                        % the last right boundary
        zc_up(end) = [];                                % Removes the last left boundary
    end
    aoi = [zc_up' zc_down'];                            % Determines the areas of interest
        
    %% Determination of PPG foot location
    b = firls(51,[0 0.9],[0 0.9*pi],'differentiator');  % Designs the differentiator FIR filter
    ppg1d = filtfilt(b,1,ppgf);                         % Obtaines the 1st derivative
    ppg2d = filtfilt(b,1,ppg1d);                        % Obtaines the 2nd derivative
    ppg3d = filtfilt(b,1,ppg2d);                        % Obtaines the 3rd derivative
    
    if v == 1
        ax(2) = subplot(2,1,2);
        yyaxis right;
        plot(ppg1d/max(ppg1d),':','color',color(7,:),'linewidth',1.2);
        plot(ppg2d/max(ppg2d),':','color',color(8,:),'linewidth',1.2);
        plot(ppg3d/max(ppg3d),'-','color',color(1,:),'linewidth',1.2);
        xlim([1 length(ppg3d)]);
        
        linkaxes(ax,'x');   
    end
    
    ons = zeros(size(aoi,1),1);                         % Initializes variable
    for i = 1:size(aoi,1)                               % Iterates through the areas of interest
        aux_ppg3d = ppg3d(aoi(i,1):aoi(i,2));           % Takes the portion of the 3rd derivative in the i-th AOI
        j = 2;                                          % Initializes counter
        zc = 0;                                         % Initializes variable
        
        while j <= length(aux_ppg3d)                    % Iterates through the portion of the signal
            if aux_ppg3d(j - 1) <= 0 && aux_ppg3d(j) > 0
                                                        % Verifies if the point corresponds to a zero-crossing
                zc = j;                                 % Updates variables
                j = length(aux_ppg3d);                  % Updates counter
            end
            j = j + 1;                                  % Updates counter
        end
        
        if zc ~= 0
            zc = zc + aoi(i,1) - 1;                     % Adds offset
            if ppg2d(zc) >= 0 && ppg1d(zc) > 0          % Verifies if the 1st and 2nd derivatives are positive at the point
                ons(i) = zc;
            else
                [~,ind] = max(ppg2d(aoi(i,1):aoi(i,2)));% Finds the location of the maximum of the second 
                                                        % derivative in the i-th AOI
                ons(i) = ind + aoi(i,1) - 1;            % Adds offset
            end
        else
            [~,ind] = max(ppg2d(aoi(i,1):aoi(i,2)));    % Finds the location of the maximum of the second 
                                                        % derivative in the i-th AOI
            ons(i) = ind + aoi(i,1) - 1;                % Adds offset
        end
    end
    
    if v == 1
        ax(1) = subplot(2,1,1);
        plot(ons,ppg(ons),'ok','MarkerFaceColor',color(1,:));
        
        ax(2) = subplot(2,1,2);
        yyaxis left;
        plot(ons,ppgf(ons),'ok','MarkerFaceColor',color(9,:));
        xlim([1 length(ppg3d)]);
        
        linkaxes(ax,'x');   
    end
    
    if v == 1
        ax(1) = subplot(2,1,1);
        miny = min(ppg);
        if miny < 0
            miny = 1.2*miny;
        else
            miny = 0.8*miny;
        end
        maxy = max(ppg);
        if maxy < 0
            maxy = 0.8*maxy;
        else
            maxy = 1.2*maxy;
        end
        yline = linspace(miny,maxy,1000);
        plot(aoi(1,1).*ones(size(yline)),yline,'--','color',color(5,:),'linewidth',1.2);
        plot(aoi(1,2).*ones(size(yline)),yline,'--','color',color(6,:),'linewidth',1.2);        
        plot(aoi(2:end,1).*ones(size(yline)),yline,'--','color',color(5,:),'linewidth',1.2);
        plot(aoi(2:end,2).*ones(size(yline)),yline,'--','color',color(6,:),'linewidth',1.2);
        legend('PPG','Onset','AOI: Left boundaries','AOI: Right boundaries','location','northeastoutside');
        
        ax(2) = subplot(2,1,2);
        yyaxis right;
        miny = min([min(ppgf) min(ppgw) min(ppgrs)]);
        if miny < 0
            miny = 1.2*miny;
        else
            miny = 0.8*miny;
        end
        maxy = max([max(ppgf) max(ppgw) max(ppgrs)]);
        if maxy < 0
            maxy = 0.8*maxy;
        else
            maxy = 1.2*maxy;
        end
        yline = linspace(miny,maxy,1000);
        plot(aoi(1,1).*ones(size(yline)),yline,'--','color',color(5,:),'linewidth',1.2);
        plot(aoi(1,2).*ones(size(yline)),yline,'--','color',color(6,:),'linewidth',1.2);
        plot(aoi(2:end,1).*ones(size(yline)),yline,'--','color',color(5,:),'linewidth',1.2);
        plot(aoi(2:end,2).*ones(size(yline)),yline,'--','color',color(6,:),'linewidth',1.2);
        legend( 'Filtered PPG','Wavelet scale','Rectified and squared','Threshold','Onset', ...
                '1D PPG','2D PPG','3D PPG','AOI: Left boundaries','AOI: Right boundaries',...
                'location','northeastoutside');
        
        linkaxes(ax,'x');  
    end
    ibis = ons;                                         % Sets the output
end

function h = waveleth(w)
%% Constructs the LPF required for the wavelet-based delineator at a sampling frequency of 250 Hz
% Inputs:   array with the frequency points, in radians, that will be used to construct the filter. 
%           Must be between 0 and 2pi
% Outputs:  array with the LPF coefficients

    aux1 = exp((1j)*w/2);
    aux2 = cos(w/2).^3;
    h = aux1.*aux2;

end

function g = waveletg(w)
%% Constructs the HPF required for the wavelet-based delineator at a sampling frequency of 250 Hz
% Inputs:   array with the frequency points, in radians, that will be used to construct the filter. 
%           Must be between 0 and 2pi
% Outputs:  array with the HPF coefficients

    aux1 = (4*(1j));
    aux2 = exp((1j)*w/2);
    aux3 = sin(w/2);
    g = (aux1*aux2).*aux3;

end

function [q, q250] = wavelet_filters(n,fs,v)
%% Creates the filters required to make the wavelet decomposition using the algorithme-a-trous
% Inputs:   number of samples of the signal that will be decomposed, n
%           sampling frequency of the signal that will be decomposed, fs
%           visualization option, v
% Outputs:  cell object with the five filters required to make the wavelet decomposition, q

    %% Generation of filters
    m = 250*n/fs;
    w = 0:2*pi/m:2*pi - 2*pi/m;
    
    q = cell(5,1);
    q{1} = waveletg(w);
    if v == 1
        figure;
        subplot(1,2,1);
        f = (0:length(q{1}) - 1)*(250/length(q{1}));
        plot(f,abs(q{1}));
        hold on;
        grid on; box off;
%         xlim([0 floor(250/2)]);
    end
    for k = 2:5
        g = waveleth(w);
        h = 1;
        while h < k - 1
            g = g.*waveleth((2^h).*w);
            h = h + 1;
        end
        g = waveletg(((2^(k - 1))*w)).*g;
        q{k} = g;
        if v == 1
            subplot(1,2,1);
            f = (0:length(q{k}) - 1)*(250/length(q{k}));
            plot(f,abs(q{k}));
            hold on;
            grid on; box off;
%           xlim([0 floor(250/2)]);
        end
    end
    q250 = q;
    
    
    %% Resampling to fs
    for k = 1:length(q)
        aux = real(ifft(q{k}));
        len = length(aux);
        xi = (0:length(aux) - 1)/250;
        xf = 0:1/fs:len/250; xf(end) = [];
        aux = spline(xi,aux,xf);
        q{k} = fft(aux); 
    end
    
    if v == 1
        subplot(1,2,2);
        hold on;
        for k = 1:length(q)
            f = (0:length(q{k}) - 1)*(fs/length(q{k}));
            plot(f,abs(q{k}));
        end
        grid on; box off;
%         xlim([0 fs/2]);
    end
end

function w = wavelet_decomposition(x,q,v)
%% Performs the Wavelet decomposition of a signal using the algorithme-a-trous
% Inputs:   signal to be decomposed, x
%           cell object containing the filters that will decompose the signal, q
%           visualization option, v
% Outputs:  cell object containing the Wavelet decomposition of the signal
%           at scales 2^1, 2^2, ..., 2^5, w

    %% Application of the Wavelet filters
    if v == 1
        figure('position',[250   90  980  680]);
        ax(1) = subplot(6,1,1);
        plot(x);
        grid on; box off;
        title('Original signal');
        linkaxes(ax,'x');
    end
    
    w = cell(size(q));
    for i = 1:length(q)
        aux = fft(x);
        aux_q = q{i};
        if size(aux,1) > size(aux,2)
            aux_q = aux_q';
        end
        aux = aux.*aux_q;
        aux = ifft(aux);
        w{i} = real(aux);
        
        if v == 1
            ax(i + 1) = subplot(6,1,i + 1);
            plot(w{i});
            grid on; box off;
            title(strcat(['Scale ' num2str(i)]));
            linkaxes(ax,'x');
        end
    end
    
end