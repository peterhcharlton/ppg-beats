function [peaks, onsets] = ppgpulses_beat_detector(sig,fs)
% PPGPULSES_BEAT_DETECTOR  PPGPULSES PPG beat detector.
%   PPGPULSES_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
%   using the 'PPG Pulses' beat detector
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
%   J. Lazaro et al., 'Pulse rate variability analysis for discrimination of sleep-apnea-related decreases in the amplitude fluctuations of pulse photoplethysmographic signal in children,' IEEE Journal of Biomedical and Health Informatics, vol. 18, no. 1, pp. 240-246, 2014. <https://doi.org/10.1109/JBHI.2013.2267096>
%   
%   # Author
%   * Jesus Lazaro - created the original code
%   * Mariano Llamedo Soria - adapted it to the ecg-kit project
%   * Peter H. Charlton - did very little, just wrote this wrapper
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Version
%   1.0
%   
%   # Source
%   This function was copied from the 'ecg-kit' repository by 'marianux' at:
%   <https://github.com/marianux/ecg-kit> , where it is licensed under the 
%   GNU License v2.0. See also: <http://marianux.github.io/ecg-kit/>
%   This particular function was available at: <https://github.com/marianux/ecg-kit/blob/master/common/ppg/PPG_pulses_detector.m>
%
%   # Licence
%      This is available under the GNU General Public License v2.0: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html

[peaks, sig_filt, peaks_filt, thres] = PPG_pulses_detector(sig, fs);

peaks = tidy_beats(peaks);

onsets = pulse_onsets_from_peaks(sig, peaks);

end

function [ peaks, sig_filt, peaks_filt, thres ] = PPG_pulses_detector( ppg, fsppg, pb, lpd_fp, lpd_fc, lpd_order, alpha, refract, taoRR, w_nA, plotflag )
% -- Source --
% This function was copied from the 'ecg-kit' repository by 'marianux' at: https://github.com/marianux/ecg-kit . See also: http://marianux.github.io/ecg-kit/ 
% This particular function was available at: https://github.com/marianux/ecg-kit/blob/master/common/ppg/PPG_pulses_detector.m
% It is available under the GNU General Public License v2.0: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
% ------------

% added by PC:
ignore_pb = 1; % skips all lines involving the 'pb' variable

%PPG_PULSES_DETECTOR     PPG signal pulses detector, based on low-pass
%                        differentiator filter.
%
% Created by Jesus Lazaro <jlazarop@unizar.es> in 2012
% adapted by Mariano Llamedo Soria to the ecg-kit project.
%--------
%   Sintax: [ peaks, sig_filt, peaks_filt, thres ] = PPG_pulses_detector( ppg, fsppg, lpd_fp, lpd_fc, lpd_order, alpha, refract, taoRR, w_nA, plotflag )
%   In:   ppg = PPG signal
%         fsppg = "ppg" sampling rate (Hz)
%         lpd_fp = last frequency of the pass-band for the
%                  low-pass-differentiator filter (Hz) [Default: 7.8]
%         lpd_fc = cut-off frequency for the low-pass-differentiator filter
%                  (Hz) [Default: 8]
%         lpd_order = filter order [Default: 3*fsppg]
%         alpha = multiplies previus amplitude of detected maximum in
%                 filtered signal for updating the threshold [Default: 0.2]
%         refract = refractary period for threshold (s) [Default: 150e-3]
%         taoRR = fraction of estimated RR where threshold reaches its
%                 minimum value (alpha*amplitude of previous SSF peak)
%                 [Default: 1]
%         w_nA = length of the window after a peak in filtered signal in
%                which is asociated peak in PPG will be searched (s)
%                [Default: 300e-3]
%         plotflag = if 1, plots a figure with PPG and SSF [Default: 0]
%
%   Out:  peaks = location of maximum of detected pulses (samples)
%         sig_filt = band-pass filtered signal
%         peaks_filt = location of peaks detected in filtered signal (samples)
%         thres = computed time varying theshold
    
    if nargin<2
        error('Not enough input arguments');
    end
    
    if nargin<3
        if ~ignore_pb
            pb = [];
        end
    end
    
    if nargin<4
        lpd_fp = 7.8;
    end
    
    if nargin<5
        lpd_fc = 8;
    end
    
    if nargin<6
        lpd_order = 3*fsppg;
    end
    
    if nargin<7
        alpha = 0.2;
    end
    
    if nargin<8
        refract = 150e-3;
    end
    
    if nargin<9
        taoRR = 1;
    end
    
    if nargin<10
         w_nA = 300e-3;
    end
    
    if nargin<11
        plotflag = 0;
    end
    
    
    ppg = ppg(:); %Ensure "ppg" is a column vector
    refract = round(refract*fsppg); %Seconds->samples
    
    %% Filtering:
    
    this_path = fileparts(mfilename('fullpath'));
    % default folder to look at
    cached_filter_filename = [this_path filesep sprintf( 'low_pass_differentiator_%dHz.mat', fsppg) ];
    if( exist(cached_filter_filename, 'file') )
        pepe = load(cached_filter_filename);
        bb = pepe.bb;
        clear pepe
    else
        % force an integer group delay
        if( rem(lpd_order,2) ~= 0 )
            lpd_order = lpd_order + 1;
        end
        d = fdesign.differentiator('n,fp,fst', lpd_order, lpd_fp*2/fsppg, lpd_fc*2/fsppg);
        % hd = design(d,'equiripple');
        hd = design(d, 'firls');
        bb = hd.Numerator*fsppg/(2*pi);
        save(cached_filter_filename, 'bb');
        clear d hd;
    end    
    
    %filter the signal
    delay = (numel(bb)-1)/2;
    sig_filt = filter(bb, 1, ppg);
    sig_filt = [sig_filt(1+delay:end); zeros(delay, 1)];
    
    
    %% Compute threshold for filtered signal:
    if ~ignore_pb
        pb.reset();
        pb.checkpoint('Compute threshold for filtered signal');
    end
    
    peaks_filt = [];
    thres_ini_w_ini = find(~isnan(sig_filt), 1, 'first');
    thres_ini_w_end = thres_ini_w_ini + round(10*fsppg);
    aux = sig_filt(thres_ini_w_ini:thres_ini_w_end);
    thres_ini = 3*mean(aux(aux>=0));
    thres = nan(size(sig_filt));
    t = 1:length(sig_filt);
    RR = round(60/80*fsppg);
    if (1+RR)<length(sig_filt)
        thres(1:1+RR) = thres_ini - (thres_ini*(1-alpha)/RR)*(t(1:RR+1)-1);
        thres(1+RR:end) = alpha*thres_ini;
    else
        thres(1:end) = thres_ini - (thres_ini*(1-alpha)/RR)*(t(1:end)-1);
    end
    
    if ~ignore_pb
        pb.Loops2Do = round(length(sig_filt) / RR);
    end
    
    kk=1;
    while true
        if ~ignore_pb
            pb.start_loop();
        end
        
        cross_u = kk-1 + find(sig_filt(kk:end)>thres(kk:end), 1, 'first'); %Next point to cross the actual threshold (down->up)
        if isempty(cross_u)
            % No more pulses -> end
            break;
        end
        
        cross_d = cross_u-1 + find(sig_filt(cross_u:end)<thres(cross_u:end), 1, 'first'); %Next point to cross the actual threshold (up->down)
        
        if isempty(cross_d)
            % No more pulses -> end
            break;
        end
        
        % Pulse detected:
        [vmax, imax] = max(sig_filt(cross_u:cross_d));
        p = cross_u-1+imax;
        peaks_filt = [peaks_filt, p];
        
        if ~ignore_pb
            pb.checkpoint([]);
        end
        
        % Update threshold
        N_RR_estimation = 3;
        N_ampli_est = 3;
        Npeaks = length(peaks_filt);
        if Npeaks>=N_RR_estimation+1;
            RR = round(median(diff(peaks_filt(end-N_RR_estimation:end))));
        elseif Npeaks>=2
            RR = round(mean(diff(peaks_filt)));
        end
        kk = min(p+refract, length(sig_filt));
        thres(p:kk) = vmax;
%         tao = 5/(taoRR*RR-refract);
%         thres(kk:end) = vmax*(1-alpha)*exp(-tao*(t(kk:end)-kk)) + vmax*alpha;
        
        if ~ignore_pb
            pb.checkpoint([]);
        end

        vfall = vmax*alpha;
        if Npeaks>=(N_ampli_est+1)
            ampli_est = median(sig_filt(peaks_filt(end-N_ampli_est:end-1)));
            if vmax>=(2*ampli_est)
                vfall = alpha*ampli_est;
                vmax = ampli_est;
            end
%             if vmax<=(0.5*ampli_est)
%                 vfall = alpha*ampli_est;
%                 vmax = ampli_est;
%             end
        end
        
        fall_end = round(taoRR*RR);
        if (kk+fall_end)<length(sig_filt)
            thres(kk:kk+fall_end) = vmax - (vmax-vfall)/fall_end*(t(kk:kk+fall_end)-kk);
            thres(kk+fall_end:end) = vfall;
        else
            thres(kk:end) = vmax - (vmax-vfall)/fall_end*(t(kk:end)-kk);
        end
        
        if ~ignore_pb
            pb.end_loop();
        end
        
    end
    
    
    %% Peak maximum search in PPG:
    if ~ignore_pb
        pb.reset();
        pb.checkpoint('Peak maximum search in PPG');
        pb.Loops2Do = length(peaks_filt);
    end
    
    peaks = nan(size(peaks_filt));
    w_nA_samples = round(fsppg*w_nA);
    for kk=1:length(peaks_filt)
        
        if ~ignore_pb
            pb.start_loop();
        end
        
%         w_int.begin = peaks_filt(kk) - w_nA_samples;
        w_int.begin = peaks_filt(kk);
        w_int.end = peaks_filt(kk) + w_nA_samples;

        if ~ignore_pb
            pb.checkpoint([]);
        end
        
        [aux, aux_t] = extract_interval(ppg, 1:length(ppg), w_int.begin, w_int.end);
        [~, pos] = max(aux);
        
        % Avoid possible detections at beggining or ending samples:
        if  ( aux_t(pos)==1 ) || ( ( aux_t(pos)==numel(ppg) ) )
            peaks(kk) = nan;
            continue;
        end
        
        if ( ppg(aux_t(pos))>=ppg(aux_t(pos)-1) ) && ( ppg(aux_t(pos))>=ppg(aux_t(pos)+1) )
            % Pulse is detected on a relative maximum -> right
            peaks(kk) = aux_t(pos);
        else
            % Pulse is not detected on a relative maximum -> wrong
            peaks(kk) = nan;
        end
        
        if ~ignore_pb
            pb.end_loop();
        end
        
    end
    peaks = peaks(~isnan(peaks));
    
    
    %% Figure:
    if plotflag==1
        figure;
        ax(1) = subplot(2,1,1); hold on;
        plot(ppg, 'b');
        plot(peaks_filt, ppg(peaks_filt), 'ro');
        plot(peaks, ppg(peaks), 'r*');
        title('PPG');
        ax(2) = subplot(2,1,2); hold on;
        plot(sig_filt, 'b');
        plot(thres, 'k');
        plot(peaks_filt, sig_filt(peaks_filt), 'ro');
        title('SSF');
        linkaxes(ax, 'x');
    end
    
end

function [ x_int, t_int, indexes ] = extract_interval( x, t, int_ini, int_end )
% -- Source --
% This function was copied from the 'ecg-kit' repository by 'marianux' at: https://github.com/marianux/ecg-kit . See also: http://marianux.github.io/ecg-kit/ 
% This particular function was available at: https://github.com/marianux/ecg-kit/blob/master/common/ppg/PPG_pulses_detector.m
% It is available under the GNU General Public License v2.0: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
% ------------

%EXTRACT_INTERVAL     Very simple function to extract an interval from a signal
%
% Created by Jesus Lazaro <jlazarop@unizar.es> in 2011
%--------
%   Sintax: [ x_int, t_int, indexes ] = extract_interval( x, t, int_ini, int_end )
%   In:   x = signal
%         t = time vector
%         int_ini = interval begin time (same units as 't')
%         int_end = interval end time (same units as 't')
%
%   Out:  x_int = interval [int_ini, int_end] of 'x'
%         t_int = interval [int_ini, int_end] of 't'
%         indexes = indexes corresponding to returned time interval

    if nargin<4
        error('Not enough input arguments');
    end
    
    indexes = find(t>=int_ini & t <=int_end);
    x_int = x(indexes);
    t_int = t(indexes);
end
