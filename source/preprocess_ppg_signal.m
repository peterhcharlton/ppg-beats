function S = preprocess_ppg_signal(S, options)
% PREPROCESS_PPG_SIGNAL  pre-processes PPG signal.
%   PREPROCESS_PPG_SIGNAL filters and downsamples a
%   photoplethysmogram (PPG) signal ready for beat detection.
%   
%   # Inputs
%   
%   * S : a structure containing three fields:
%    - v : a vector of PPG values
%    - fs : the sampling frequency of the PPG in Hz
%    - no_signal : (optional) a logical indicating whether (1) or not (0) each PPG sample contains a physiological signal
%    
%   * options - A structure containing the settings used for the processing:
%    - options.do_downsample   - (1 or 0, default value of 0) A logical indicating whether or not to downsample the PPG signals prior to analysis.
%    - options.downsample_freq - the frequency at which to downsample PPG signals (Hz). Only used if downsampling is enabled by options.do_downsample.
%    - options.redo_analysis   - A logical indicating whether or not to overwrite existing analysis files in order to redo the analysis
%    - options.filtering.elim_high_freqs - ...TBC...
%    - options.filtering.elim_low_freqs  - ...TBC...
%   
%   # Outputs
%   * S : containing the filtered and downsampled PPG signal
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Author
%   Peter H. Charlton, University of Cambridge, February 2022.
%   
%   # Version
%   1.0
%   
%   # License - GPL-3.0
%      Copyright (c) 2022 Peter H. Charlton
%      This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%      This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%      You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

% Downsample 'no_signal'
if options.do_downsample && sum(strcmp(fieldnames(S), 'no_signal'))
    if rem(S.fs, options.downsample_freq) == 0
        % Downsample
        ds_factor = S.fs/options.downsample_freq;
        S.no_signal = downsample(S.no_signal, ds_factor);
    else
        % Resample
        %S.no_signal = resample(double(S.no_signal), options.downsample_freq, S.fs);
        t = [0:length(S.no_signal)-1]./S.fs;
        S.no_signal = interp1(t, double(S.no_signal), 0:(1/options.downsample_freq):t(end), 'previous');
        S.no_signal = logical(S.no_signal);
    end
end

% Downsample PPG signal
if options.do_downsample
    S = do_downsample(S, options.downsample_freq);
end

% Filter PPG signal
sigs = filter_signal(S, options);
S.v = sigs.filt;
S.fs = sigs.fs;
clear options

end

function S = do_downsample(S, downsample_freq)

%% Downsample

% Work out whether to downsample or interpolate
if rem(S.fs, downsample_freq) == 0
    % Downsample
    ds_factor = S.fs/downsample_freq;
    S.v = downsample(S.v, ds_factor);
    S.fs = downsample_freq;
else
    use_resampling = 1;
    if use_resampling
        % Resample
        S.v(isnan(S.v)) = median(S.v(~isnan(S.v)));
        S.v = resample(S.v, downsample_freq, S.fs);
    else
        % Interpolate
        t = 1:length(S.v);
        t_new = 1:(S.fs/downsample_freq):t(end);
        S.v = interp1(t, S.v, t_new);
    end
    S.fs = downsample_freq;
end

end

function sigs = filter_signal(S, options)

%% Make output variable
sigs.fs = S.fs;
sigs.orig = S.v;

%% Filter to remove high frequencies
filt_characteristics = options.filtering.elim_high_freqs;
s_filt = elim_vhfs3(S, filt_characteristics);

%% Filter to remove low frequencies
filt_characteristics = options.filtering.elim_low_freqs;
s_filt = elim_vlfs(s_filt, filt_characteristics);

%% Output
sigs.filt = s_filt.v;
sigs.curr = sigs.filt;

end

function s_filt = elim_vhfs3(s, filt_characteristics)
%% Filter signal to remove VHFs
% Adapted from RRest

s_filt.fs = s.fs;

%% Eliminate nans
s.v(isnan(s.v)) = mean(s.v(~isnan(s.v)));

%% Check to see if sampling freq is at least twice the freq of interest
if (filt_characteristics.Fpass/(s.fs/2)) >= 1
    % then the fs is too low to perform this filtering
    s_filt.v = s.v;
    return
end

%% Create filter
% parameters for the low-pass filter to be used
AMfilter = create_lpf(filt_characteristics, s);

%% Re-make filter if it requires too many samples for this signal
% check to see if it requires too many samples
req_sig_length = 3*(length(AMfilter.numerator)-1);
no_attempts = 0;
while no_attempts<4 && length(s.v)<=req_sig_length
    % change Fpass (i.e. the high frequency end of the filter band)
    filt_characteristics.Fpass = filt_characteristics.Fpass+0.75*(filt_characteristics.Fpass-filt_characteristics.Fstop);
    % re-make filter
    AMfilter = create_lpf(filt_characteristics, s);
    % update criterion
    req_sig_length = 3*(length(AMfilter.numerator)-1);
    % increment number of attempts
    no_attempts = no_attempts+1;
end
if length(s.v)<=req_sig_length
    fprintf('\n - Couldn''t perform high frequency filtering')
end

%% Check frequency response
% Gives a -3 dB cutoff at cutoff_freq Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 0.0512;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(s.fs/2);

%% Remove VHFs
if length(s.v)>req_sig_length
    s_filt.v = filtfilt(AMfilter.numerator, 1, s.v);
else
    s_filt.v = s.v;
end

end

function AMfilter = create_lpf(filt_characteristics, s)

[N,Wn,BETA,TYPE] = kaiserord([filt_characteristics.Fstop filt_characteristics.Fpass]/(s.fs/2), [1 0], [filt_characteristics.Dstop filt_characteristics.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), 'scale');
AMfilter = dfilt.dffir(b);

end

function s_filt = elim_vlfs(s, filt_characteristics)
%% Filter pre-processed signal to remove frequencies below resp
% Adapted from RRest

%% Eliminate nans
s.v(isnan(s.v)) = mean(s.v(~isnan(s.v)));

%% Make filter
flag  = 'scale';
[N,Wn,BETA,TYPE] = kaiserord([filt_characteristics.Fstop filt_characteristics.Fpass]/(s.fs/2), [1 0], [filt_characteristics.Dstop filt_characteristics.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% % Gives a -3 dB cutoff at ? Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 4.435e-3;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(s.fs/2);

if length(s.v) > (length(AMfilter.numerator)-1)*3
    % - Tukey Window to avoid edge effects
    win_edge_durn = 0.5; % in secs
    prop_win = win_edge_durn*2/((length(s.v)-1)/s.fs);
    tw = tukeywin(length(s.v),prop_win); 
    s_filt.v = filtfilt(AMfilter.numerator, 1, s.v.*tw);
    s_filt.v = s.v-s_filt.v;
else
    s_filt.v = s.v;
end
s_filt.fs = s.fs;
end