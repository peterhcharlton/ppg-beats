function make_plot_of_ppg_beat_detection_challenges

% make subplot for each challenge
challenges = {'sitting', 'working', 'cycling', 'walking', 'lunch_break', 'car_driving', 'stair_climbing', 'table_soccer'}; %, 'af'};
for challenge_no = 1 : length(challenges)
    curr_challenge = challenges{challenge_no};
    
    %% Load raw data

    % setup paths for this challenge
    subj_no = 2;
    switch curr_challenge
        case 'sitting'
            raw_data_path = '/Users/petercharlton/Documents/Data/ppg_dalia/conv_data/ppg_dalia_sitting_data.mat';
            
        case 'working'
            raw_data_path = '/Users/petercharlton/Documents/Data/ppg_dalia/conv_data/ppg_dalia_working_data.mat';
            
        case 'cycling'
            raw_data_path = '/Users/petercharlton/Documents/Data/ppg_dalia/conv_data/ppg_dalia_cycling_data.mat';
            
        case 'walking'
            raw_data_path = '/Users/petercharlton/Documents/Data/ppg_dalia/conv_data/ppg_dalia_walking_data.mat';
            
        case 'lunch_break'
            raw_data_path = '/Users/petercharlton/Documents/Data/ppg_dalia/conv_data/ppg_dalia_lunch_break_data.mat';
            
        case 'car_driving'
            raw_data_path = '/Users/petercharlton/Documents/Data/ppg_dalia/conv_data/ppg_dalia_car_driving_data.mat';
            
        case 'stair_climbing'
            raw_data_path = '/Users/petercharlton/Documents/Data/ppg_dalia/conv_data/ppg_dalia_stair_climbing_data.mat';
            
        case 'table_soccer'
            raw_data_path = '/Users/petercharlton/Documents/Data/ppg_dalia/conv_data/ppg_dalia_table_soccer_data.mat';
            
    end
    
    % load raw data for this challenge
    load(raw_data_path);
    
    % Obtain beats for this challenge
    load_beats = 0;
    if load_beats
        % Load PPG beats 
        % - Identify path of detected PPG beats for this subject
        subj_no_txt = num2str(subj_no, '%04.f');
        temp = strfind(raw_data_path, '/');
        temp2 = strfind(raw_data_path, '_');
        ppg_beats_path = [raw_data_path(1:temp(end)), 'proc_data_', raw_data_path(temp(end)+1:temp2(end)-1), filesep, subj_no_txt, '_ppg_beats'];
        load(ppg_beats_path);
        ecg_beats_aligned_path = [raw_data_path(1:temp(end)), 'proc_data_', raw_data_path(temp(end)+1:temp2(end)-1), filesep, subj_no_txt, '_ecg_beats_aligned'];
        load(ecg_beats_aligned_path);
    else
        % Identify beats
        [ppg_beats_inds.MSPTD, ~, ~] = detect_ppg_beats(data(subj_no).ppg, 'MSPTD');     % detect beats in PPG
        [ppg_beats_inds.qppgfast, ~, ~] = detect_ppg_beats(data(subj_no).ppg, 'qppgfast');     % detect beats in PPG
        [ecg_beats_inds, qual] = detect_ecg_beats(data(subj_no).ecg.v, data(subj_no).ecg.fs);           % detect beats in ECG and assess quality of beat detection
        ecg_exc_log = ~qual;
        options = struct;
        [ecg_beats_a_inds.MSPTD, ecg_exc_a_log, lag_ecg_samps] = align_ppg_ecg_beats(ppg_beats_inds.MSPTD, ecg_beats_inds, data(subj_no).ppg.fs, data(subj_no).ecg.fs, options, ecg_exc_log);
        
    end
    
    % extract required ECG and PPG data for this challenge
    rel_data(challenge_no).ecg = data(subj_no).ecg;
    rel_data(challenge_no).ecg.t = [0:length(rel_data(challenge_no).ecg.v)-1]./rel_data(challenge_no).ecg.fs;
    rel_data(challenge_no).start_time = find(ecg_exc_a_log==0,1)/data(subj_no).ecg.fs;
    rel_data(challenge_no).ppg = data(subj_no).ppg; clear data
    rel_data(challenge_no).ppg.t = [0:length(rel_data(challenge_no).ppg.v)-1]./rel_data(challenge_no).ppg.fs;
    rel_data(challenge_no).ppg_beats_MSPTD = ppg_beats_inds.MSPTD;
    rel_data(challenge_no).ppg_beats_qppg = ppg_beats_inds.qppgfast;
    rel_data(challenge_no).ecg_beats = ecg_beats_a_inds.MSPTD;
    
end


% - setup plot
close all
figure('Position', [20,20,1162,930])
ftsize = 17;
lwidth = 1.5;

% - make a subplot for each challenge
for challenge_no = 1 : length(challenges)
    curr_challenge = challenges{challenge_no};
    
    % setup title
    switch curr_challenge
        case 'sitting'
            txt = '(a) Sitting';
        case 'working'
            txt = '(b) Working';
        case 'cycling'
            txt = '(c) Cycling';
        case 'walking'
            txt = '(d) Walking';
        case 'lunch_break'
            txt = '(e) Lunch break';
        case 'car_driving'
            txt = '(f) Car driving';
        case 'stair_climbing'
            txt = '(g) Stair climbing';
        case 'table_soccer'
            txt = '(h) Table soccer';
    end
    
    %% Process data for plot
    
    
    % - select a 10-second window
    durn = 10; % duration in seconds
    end_time = rel_data(challenge_no).start_time + durn;
    rel_els = find(rel_data(challenge_no).ppg.t >= rel_data(challenge_no).start_time & rel_data(challenge_no).ppg.t <= end_time);
    rel_ppg.fs = rel_data(challenge_no).ppg.fs;
    rel_ppg.v = rel_data(challenge_no).ppg.v(rel_els);
    rel_ppg.t = [0:length(rel_ppg.v)-1]./rel_ppg.fs;
    rel_ppg.MSPTD = rel_data(challenge_no).ppg_beats_MSPTD(rel_data(challenge_no).ppg_beats_MSPTD>= rel_els(1) & rel_data(challenge_no).ppg_beats_MSPTD<= rel_els(end)) - rel_els(1)+1;
    rel_ppg.qppg = rel_data(challenge_no).ppg_beats_qppg(rel_data(challenge_no).ppg_beats_qppg>= rel_els(1) & rel_data(challenge_no).ppg_beats_qppg<= rel_els(end)) - rel_els(1)+1;
    rel_ecg = rel_data(challenge_no).ecg;
    rel_els = find(rel_ecg.t >= rel_data(challenge_no).start_time & rel_ecg.t <= end_time);
    rel_ecg.beats = rel_data(challenge_no).ecg_beats(rel_data(challenge_no).ecg_beats>= rel_els(1) & rel_data(challenge_no).ecg_beats<= rel_els(end)) - rel_els(1)+1;
    
    % - normalise PPG signal to lie between 0 and 1
    rel_ppg.v = (rel_ppg.v-min(rel_ppg.v))./range(rel_ppg.v);

    %% Make plot
    %subplot(length(challenges),1, challenge_no)
    if length(challenges) == 4
        subplot('Position', [0.08,0.06+0.24*(length(challenges)-challenge_no), 0.89, 0.18])
    else
        curr_y_number = 4 - challenge_no;
        if curr_y_number < 0
            curr_y_number = curr_y_number+4;
        end
        if challenge_no < 5
            subplot('Position', [0.08,0.06+0.24*(curr_y_number), 0.40, 0.18])
        else
            subplot('Position', [0.59,0.06+0.24*(curr_y_number), 0.40, 0.18])
        end
    end
    
    % - plot ECG beats
    for beat_no = 1 : length(rel_ecg.beats)
        h = plot(ones(2,1)*rel_ecg.t(rel_ecg.beats(beat_no)), [-0.1,0.23], '--', 'Color', 0.2*ones(1,3), 'LineWidth', lwidth);
        hold on
    end
    
    % - plot PPG signal
    plot(rel_ppg.t, rel_ppg.v, 'b', 'LineWidth', lwidth), hold on

    % - plot PPG beats
    beat_detectors = {'MSPTD', 'qppg'};
    mk_styles = {'+k','or'};
    mk_sizes = [16,12];
    for beat_detector_no = 1:length(beat_detectors)
        curr_beat_detector = beat_detectors{beat_detector_no};
        eval(['rel_ppg_beats_inds = rel_ppg.' curr_beat_detector ';']);
        h(beat_detector_no+1) = plot(rel_ppg.t(rel_ppg_beats_inds), rel_ppg.v(rel_ppg_beats_inds), mk_styles{beat_detector_no}, 'LineWidth', lwidth, 'MarkerSize', mk_sizes(beat_detector_no));
    end
    
    % - Label axes
    ylabel('PPG (au)', 'FontSize', ftsize)
    if curr_y_number ==0
        xlabel('Time (s)', 'FontSize', ftsize)
    end
    
    % - tidy up plot
    xlim([0 durn])
    ylim([-0.1, 1.1])
    xticks = 0:2:durn;
    set(gca, 'FontSize', ftsize, 'XTick', xticks, 'YTick', [], 'XGrid', 'on')
    box off
    
    % - add title
    if length(challenges) == 4
        x_val = .44;
    else
        if challenge_no < 5
            x_val = 0.2;
        elseif challenge_no >=5
            x_val = 0.8;
        end
    end
%     if curr_y_number ~=1
        dim = [x_val 0.17+0.24*(curr_y_number) .1 .1];
%     else
%         dim = [x_val 0.14+0.24*(curr_y_number) .1 .1];
%     end
    annotation('textbox',dim,'String',txt,'FitBoxToText','on', 'FontSize', ftsize, 'LineStyle', 'none');

    % legend
    if challenge_no == 1
        legend(h, {'ECG beats', 'MSPTD beat detector', 'qppg beat detector'}, 'Position', [0.525,0.98,0.01,0.01], 'Orientation', 'horizontal', 'FontSize', ftsize)
    end

end

% save figure
filepath = '/Users/petercharlton/Google Drive/Work/Images/PPG beat detection/ppg_beat_detection_challenges_new';
save_fig(filepath)

end

function save_fig(filepath)

print(gcf, filepath, '-depsc')
fid = fopen([filepath, '.txt'], 'w');
fprintf(fid, ['Created using ' mfilename, ', ', date]);
fclose(fid);

end