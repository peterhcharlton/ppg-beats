function create_ppg_beat_detection_fig(data)

close all

% load data
filepath = '/Users/petercharlton/Documents/Data/mimiciii_ppg_neonate_adult_beat_detection/mimic_train_a_data.mat';
if nargin==0
    load(filepath);
end

% create assessment framework plot
create_assessment_framework_plot(data);

% select signal of interest
subj_no = 3;
durn = 10;
buffer = 5;
S_orig = data(subj_no).ppg;
rel_els = 1:(1+durn*S_orig.fs+2*buffer*S_orig.fs);
rel_els2 = rel_els + 1000;
rel_els1 = rel_els + 49800;

% make figure
figure('Position', [20,20,500,400])

% cycle through each sample
for sample_no = 1 :2
    
    eval(['curr_rel_els = rel_els' num2str(sample_no) ';']);
    S.v = S_orig.v(curr_rel_els);
    S.fs = S_orig.fs;
    
    %% - filter signal
    
    % - Fiducial Points: eliminating high freqs
    uParams.analysis.filtering.elim_high_freqs.Fpass = 10;   % in HZ
    uParams.analysis.filtering.elim_high_freqs.Fstop = 6.93;   % in HZ   % 6.93 - 10 gives -3dB cutoff of 8.00 Hz
    uParams.analysis.filtering.elim_high_freqs.Dpass = 0.05;
    uParams.analysis.filtering.elim_high_freqs.Dstop = 0.01;
    
    % - Fiducial Points: eliminating low freqs
    uParams.analysis.filtering.elim_low_freqs.Fpass = 1.02;  % in Hz
    uParams.analysis.filtering.elim_low_freqs.Fstop = 0.6;  % in Hz      % 0.6 - 1.02 gives -3dB cutoff of 40 bpm (i.e. 0.67 Hz)
    uParams.analysis.filtering.elim_low_freqs.Dpass = 0.05;
    uParams.analysis.filtering.elim_low_freqs.Dstop = 0.01;
    
    % - setup PulseAnalyse options for filtering
    options.elim_high_freqs = uParams.analysis.filtering.elim_high_freqs;
    options.elim_low_freqs = uParams.analysis.filtering.elim_low_freqs;
    options.calc_pw_inds = 0;
    options.do_plot = 0;
    options.do_beats = 0;
    options.do_quality = 0;
    options.close_figures = 0;
    
    % Downsample and filter PPG signal
    [~, ~, ~, sigs] = PulseAnalyse(S, options);
    S.v = sigs.filt;
    S.fs = sigs.fs;
    clear options
    
    % specify beat detectors
    beat_detectors = {'MSPTD', 'Pulses'};
    
    % detect beats
    for beat_detector_no = 1 : length(beat_detectors)
        curr_detector = beat_detectors{beat_detector_no};
        
        %% - beat detection
        
        % - setup PulseAnalyse options for beat detection
        options.calc_pw_inds = 0;
        options.do_plot = 0;
        options.do_beat_filter = 0;
        options.do_filter = 0;
        options.do_quality = 0;
        options.close_figures = 0;
        options.beat_detector = curr_detector;
        
        % - detect beats in this subject's PPG signals using this peak detector
        [~, ~, pulses, ~] = PulseAnalyse(S, options);
        peaks{beat_detector_no,1} = pulses.peaks;
        clear options pulses
        
    end
    
    %% create plot
    
    % settings
    ftsize = 16;
    markers = {'+', 'o'};
    marker_sizes = [12,8];
    colors = {'k', 'r'};
    
    % make time vector
    t = [0:length(S.v)-1]/S.fs; t = t-5;
    
    % plot PPG
    if sample_no == 1
        subplot('Position', [0.08,0.56, 0.89, 0.43])
    else
        subplot('Position', [0.08,0.13, 0.89, 0.35])
    end
    plot(t, S.v, 'Color', 0.2*ones(1,3), 'LineWidth', 2)
    hold on
    
    % plot detected beats
    for detector_no = 1 : length(peaks)
        curr_peaks = peaks{detector_no};
        curr_marker = markers{detector_no};
        curr_marker_size = marker_sizes(detector_no);
        curr_color = colors{detector_no};
        h(detector_no) = plot(t(curr_peaks), S.v(curr_peaks), curr_marker, 'LineWidth', 1.5, 'Color', curr_color, 'MarkerSize', curr_marker_size);
    end
    
    % tidy up
    if sample_no == 2
        xlabel('Time (s)', 'FontSize', ftsize)
    end
    ylabel('PPG (au)', 'FontSize', ftsize)
    set(gca, 'YTick', [], 'FontSize', ftsize, 'XGrid', 'on')
    box off
    xlim([0 10])
    
    % legend
    if sample_no == 1
        legend(h, {'Beat detector 1', 'Beat detector 2'}, 'Location', 'northoutside', 'Orientation', 'horizontal')
    end
    
end

% annotate
dim = [.0 .84 .1 .1];
str = '(a)';
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none', 'FontSize', ftsize);
dim = [.0 .4 .1 .1];
str = '(b)';
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineStyle', 'none', 'FontSize', ftsize);

% save figure
filepath = '/Users/petercharlton/Google Drive/Work/Images/PPG beat detection/beat_detection';
save_fig(filepath)


end

function save_fig(filepath)

print(gcf, filepath, '-depsc')
fid = fopen([filepath, '.txt'], 'w');
fprintf(fid, ['Created using ' mfilename, ', ', date]);
fclose(fid);

end

function create_assessment_framework_plot(data)

% make figure
figure('Position', [20,20,500,400])
subplot('Position', [0.01,0.01,0.98,0.98])
lwidth = 2;

% set up
no_sigs = 10;
no_els = data(1).ppg.fs*5;
scale_vector = linspace(1, 0.2, no_els);
offset_vector = linspace(0,2, no_els);
subj_nos = [19,20,22,25,26, 29:32,34];
start_el = 125*5*60;
rel_els = start_el:start_el+no_els-1;

% extract data for each signal
for s = 1 : no_sigs
    rel_data(s).v = data(subj_nos(s)).ppg.v(rel_els);
    rel_data(s).v = (rel_data(s).v-min(rel_data(s).v))./range(rel_data(s).v);
    rel_data(s).v = rel_data(s).v(:) .* scale_vector(:) + offset_vector(:);
end

% plot data
for s = 1 : no_sigs
    plot(s+rel_data(s).v, 'b', 'LineWidth', lwidth);
    hold on
end

% tidy up
box off
xlim([0 no_els])
ylim([1 no_sigs+2.2])
set(gca, 'XTick', [], 'YTick', [], 'visible', 'off')

% save figure
filepath = '/Users/petercharlton/Google Drive/Work/Images/PPG beat detection/assessment_framework';
save_fig(filepath)
close all

end