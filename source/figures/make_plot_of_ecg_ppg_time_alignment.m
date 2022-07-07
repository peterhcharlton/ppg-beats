function make_plot_of_ecg_ppg_time_alignment

close all

% Load raw data
load('/Users/petercharlton/Documents/Data/mimiciii_ppg_neonate_adult_beat_detection/mimic_test_all_data.mat');

% specify subject
subj_no = 2; % 2 is best so far
subj_no_txt = num2str(subj_no, '%04.f');

% Load beats data
load(['/Users/petercharlton/Documents/Data/mimiciii_ppg_neonate_adult_beat_detection/proc_data_mimic_test_all/' subj_no_txt '_ecg_beats.mat'])
load(['/Users/petercharlton/Documents/Data/mimiciii_ppg_neonate_adult_beat_detection/proc_data_mimic_test_all/' subj_no_txt '_ppg_beats.mat'])

% Extract PPG and ECG signals
ecg = data(subj_no).ekg;
ppg = data(subj_no).ppg;

% Detect beats in two signals
beat_detector = 'qppgfast';     
[peaks, onsets, mid_amps] = detect_ppg_beats(ppg, beat_detector);     % detect beats in PPG
ppg_beats_inds = mid_amps;
[ecg_beats_inds, qual] = detect_ecg_beats(ecg.v, ecg.fs);           % detect beats in ECG and assess quality of beat detections
ecg_exc_log = ~qual;

% Find out time-alignment
options = struct;
[ecg_beats_a_inds, ecg_exc_a_log, lag_ecg_samps] = align_ppg_ecg_beats(ppg_beats_inds, ecg_beats_inds, ppg.fs, ecg.fs, options, ecg_exc_log);
lag_t = lag_ecg_samps/ecg.fs;
                
%% Prepare data for plot

% - calcualate time vectors
ppg.t = [0:length(ppg.v)-1]./ppg.fs;
ecg.t = [0:length(ecg.v)-1]./ecg.fs;

% - select a 5-second window
durn = 5; % duration in seconds
rel_el = find(ecg_exc_log==0,1);
rel_el = find(ecg_exc_log==0); rel_el = rel_el(42950); rel_el = 64200;
start_time = ecg.t(rel_el);
end_time = start_time + durn;
ecg_buffer = 1;
ecg_start_time = start_time - ecg_buffer;
ecg_end_time = end_time + ecg_buffer;
rel_ppg.v = ppg.v(ppg.t >= start_time & ppg.t <= end_time);
rel_ecg.v = ecg.v(ecg.t >= ecg_start_time & ecg.t <= ecg_end_time);

% - normalise signals to lie between 0 and 1 (for PPG) and 1.75 and 2.75 (for ECG)
gap = 0.5;
ppg.v = (ppg.v-min(rel_ppg.v))./range(rel_ppg.v);
ecg.v = (1+gap) + (ecg.v-min(rel_ecg.v))./range(rel_ecg.v);

%% Make plot of time alignment

do_time_alignment = 1;

if do_time_alignment
    
    % - setup
    figure('Position', [20,20,800,400])
    subplot('Position', [0.08,0.19,0.91,0.8])
    ftsize = 24;
    lwidth = 1.5;
    mksize = 12;
    
    % - plot signals
    rel_ppg_els = ppg.t >= start_time & ppg.t <= end_time;
    plot(ppg.t(rel_ppg_els), ppg.v(rel_ppg_els), 'b', 'LineWidth', lwidth), hold on
    rel_ecg_els = ecg.t >= ecg_start_time & ecg.t <= ecg_end_time;
    plot(ecg.t(rel_ecg_els), ecg.v(rel_ecg_els), 'b', 'LineWidth', lwidth)
    
    % - plot beats
    rel_ppg_beats_inds = ppg_beats_inds(ppg.t(ppg_beats_inds)>=start_time & ppg.t(ppg_beats_inds)<=end_time);
    plot(ppg.t(rel_ppg_beats_inds), ppg.v(rel_ppg_beats_inds), 'or', 'LineWidth', lwidth, 'MarkerSize', mksize)
    rel_ecg_beats_inds = ecg_beats_inds(ecg.t(ecg_beats_inds)>=(start_time-lag_t) & ecg.t(ecg_beats_inds)<=(end_time-lag_t));
    plot(ecg.t(rel_ecg_beats_inds), ecg.v(rel_ecg_beats_inds), 'ok', 'LineWidth', lwidth, 'MarkerSize', mksize)
    
    % Label axes
    xlabel('Time (s)', 'FontSize', ftsize)
    %ylabel('Signals', 'FontSize', ftsize)
    
    % tidy up plot
    xlim([ecg_start_time ecg_end_time+1])
    ylim([-0.25, 2.25+gap])
    xticks = ecg_start_time:ecg_end_time+1;
    xtick_labels = strsplit(num2str(xticks-ecg_start_time-ecg_buffer));
    set(gca, 'FontSize', ftsize, 'XTick', xticks, 'YTick', [0.5, 1.5+gap], 'YTickLabel', {'PPG', 'ECG'}, 'XTickLabel', xtick_labels, 'XGrid', 'on')
    box off
    
    
    % make dashed lines at locations of beats
    sigs = {'ecg', 'ppg'};
    for sig_no = 1 : length(sigs)
        curr_sig = sigs{sig_no};
        eval(['curr_beats = rel_' curr_sig '_beats_inds;']);
        eval(['curr_t = ' curr_sig '.t;']);
        switch curr_sig
            case 'ppg'
                ylims = [0, 1+0.5*gap];
                sig_color = 'r';
                curr_beats = curr_beats(curr_t(curr_beats) >= start_time & curr_t(curr_beats) <= end_time);
            case 'ecg'
                ylims = [1+0.5*gap, 2+gap];
                sig_color = 0.2*ones(1,3);
                curr_beats = curr_beats(curr_t(curr_beats) >= ecg_start_time & curr_t(curr_beats) <= ecg_end_time);
        end
        
        for beat_no = 1 : length(curr_beats)
            % plot dashed lines
            line_t = curr_t(curr_beats(beat_no));
            plot(ones(1,2)*line_t, ylims, '--', 'Color', sig_color, 'LineWidth', lwidth);
            
            % plot arrows
            if strcmp(curr_sig, 'ecg')
                plot([line_t, line_t+lag_t], ones(2,1)*ylims(1), 'k', 'LineWidth', lwidth)
                plot([line_t+lag_t-0.1, line_t+lag_t], [ylims(1)+0.05, ylims(1)], 'k', 'LineWidth', lwidth)
                plot([line_t+lag_t-0.1, line_t+lag_t], [ylims(1)-0.05, ylims(1)], 'k', 'LineWidth', lwidth)
                
                % annotate arrow
                if beat_no == 1
                    dim = [.10 .56 .1 .1];
                    str = 'Lag';
                    annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle', 'none', 'FontSize', ftsize - 3);
                    dim = [.09 .49 .1 .1];
                    str = [num2str(lag_t, '%.2f') 's'];
                    annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle', 'none', 'FontSize', ftsize - 3);
                end
            end
            
        end
    end
    
    % save figure
    filepath = '/Users/petercharlton/Google Drive/Work/Publications/In Preparation/2021 Assess PPG beat detectors/figures/ecg_ppg_time_alignment_new';
    save_fig(filepath)
    close all
    
end

%% Make plot of detecting correct beats

% - setup
figure('Position', [20,20,800,400])
subplot('Position', [0.08,0.19,0.91,0.8])
ftsize = 24;
lwidth = 1.5;
mksize = 12;

% - plot signals
rel_ppg_els = ppg.t >= start_time & ppg.t <= end_time;
plot(ppg.t(rel_ppg_els), ppg.v(rel_ppg_els), 'b', 'LineWidth', lwidth), hold on
rel_ecg_els = ecg.t >= ecg_start_time & ecg.t <= ecg_end_time;
plot(ecg.t(rel_ecg_els)+lag_t, ecg.v(rel_ecg_els), 'b', 'LineWidth', lwidth)

% - plot beats
rel_ppg_beats_inds = ppg_beats_inds(ppg.t(ppg_beats_inds)>=start_time & ppg.t(ppg_beats_inds)<=end_time);
plot(ppg.t(rel_ppg_beats_inds), ppg.v(rel_ppg_beats_inds), 'or', 'LineWidth', lwidth, 'MarkerSize', mksize)
rel_ecg_beats_inds = ecg_beats_inds(ecg.t(ecg_beats_inds)>=(start_time-lag_t) & ecg.t(ecg_beats_inds)<=(end_time-lag_t));
plot(ecg.t(rel_ecg_beats_inds)+lag_t, ecg.v(rel_ecg_beats_inds), 'ok', 'LineWidth', lwidth, 'MarkerSize', mksize)

% Label axes
xlabel('Time (s)', 'FontSize', ftsize)
%ylabel('Signals', 'FontSize', ftsize)

% tidy up plot
xlim([ecg_start_time ecg_end_time+1])
ylim([-0.25, 2.25+gap])
xticks = ecg_start_time:ecg_end_time+1;
xtick_labels = strsplit(num2str(xticks-ecg_start_time-ecg_buffer));
set(gca, 'FontSize', ftsize, 'XTick', xticks, 'YTick', [0.5, 1.5+gap], 'YTickLabel', {'PPG', 'ECG'}, 'XTickLabel', xtick_labels, 'XGrid', 'on')
dim = [.0 .85 .077 .1];
str = {'time-', 'aligned'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle', 'none', 'FontSize', ftsize - 6, 'HorizontalAlignment', 'center');
box off

% make dashed lines at locations of beats
sigs = {'ecg', 'ppg'};
for sig_no = 1 : length(sigs)
    curr_sig = sigs{sig_no};
    eval(['curr_beats = rel_' curr_sig '_beats_inds;']);
    eval(['curr_t = ' curr_sig '.t;']);
    if strcmp(curr_sig, 'ecg')
        curr_t = curr_t+lag_t;
    end
    switch curr_sig
        case 'ppg'
            ylims = [0, 1+0.5*gap];
            sig_color = 'r';
            curr_beats = curr_beats(curr_t(curr_beats) >= start_time & curr_t(curr_beats) <= end_time);
        case 'ecg'
            ylims = [1+0.5*gap, 2+gap];
            sig_color = 0.2*ones(1,3);
            curr_beats = curr_beats(curr_t(curr_beats) >= ecg_start_time+lag_t & curr_t(curr_beats) <= ecg_end_time+lag_t);
    end
        
    for beat_no = 1 : length(curr_beats)
        % plot dashed lines
        line_t = curr_t(curr_beats(beat_no));
        plot(ones(1,2)*line_t, ylims, '--', 'Color', sig_color, 'LineWidth', lwidth);
        % annotate
        
        
        % plot arrows
        if strcmp(curr_sig, 'ecg')
            tol = 0.15;
            plot([line_t-tol, line_t+tol], ones(2,1)*ylims(1), 'k', 'LineWidth', lwidth)
            len = 0.075;
            plot([line_t+tol-len, line_t+tol], [ylims(1)+0.05, ylims(1)], 'k', 'LineWidth', lwidth)
            plot([line_t+tol-len, line_t+tol], [ylims(1)-0.05, ylims(1)], 'k', 'LineWidth', lwidth)
            plot([line_t-tol+len, line_t-tol], [ylims(1)+0.05, ylims(1)], 'k', 'LineWidth', lwidth)
            plot([line_t-tol+len, line_t-tol], [ylims(1)-0.05, ylims(1)], 'k', 'LineWidth', lwidth)
            
            % annotate arrow
            if beat_no == 1
                dim = [.105 .56 .1 .1];
                str = 'Tol.';
                annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle', 'none', 'FontSize', ftsize - 4);
                dim = [.085 .51 .1 .1];
                str = [num2str(tol*1000) 'ms'];
                annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle', 'none', 'FontSize', ftsize - 4);
            end
        end
        
    end
end

% save figure
filepath = '/Users/petercharlton/Google Drive/Work/Publications/In Preparation/2021 Assess PPG beat detectors/figures/ppg_beat_classification_new';
save_fig(filepath)

end

function save_fig(filepath)

print(gcf, filepath, '-depsc')
fid = fopen([filepath, '.txt'], 'w');
fprintf(fid, ['Created using ' mfilename, ', ', date]);
fclose(fid);

end