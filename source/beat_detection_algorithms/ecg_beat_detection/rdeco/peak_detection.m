function [R_peak,RR_int,parameters,check] = peak_detection(parameters,signal,fs,parameters_check)
% Output:
% R_peak     - Location of the R-peaks in samples
% RR_int     - Intervals between the R-peaks in ms
% check      - Check if the method was run correctly
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

% Set check
check = 1;

% Get the size of the signal
[original_signal_length,nr_ch] = size(signal);

% Pre-allocate the R-peak and RR-interval variable
R_peak = cell(1,nr_ch);
RR_int = cell(1,nr_ch);

if parameters_check
    % Get the different parameters
    parameters = Parameterdialog(parameters,signal,fs);
end

% Check validity of the parameter dialog
if ~isempty(parameters)      
    env = round(fs*parameters{1}/1000);
    avgHR = parameters{2};
    postproc = parameters{3};
    ect = parameters{4};
    inverted = parameters{5};
    
    % Check if the signal is inverted, if so, take action
    if inverted
        signal = -signal;
    end

    % Set the segment size to one minute and compute the amount of minutes in
    % the signal
    segm_length = 60*fs; 
    nr_segm_original = floor(original_signal_length/segm_length); 

    % Add zeros if the signal does not contain a round number of minutes
    if nr_segm_original ~= 0
        signal = [signal; zeros(segm_length-(original_signal_length-nr_segm_original*segm_length),nr_ch)];
        
        % Get the new length and amount of segments
        [new_signal_length,~] = size(signal);
        nr_segm_new = floor(new_signal_length/segm_length);
    else
        % Get the new length and amount of segments
        [new_signal_length,~] = size(signal);
        nr_segm_new = 1;
    end
    
    % Pre-allocate the R-peaks per channel
    R_peak_ch = cell(nr_segm_new);

    % Loop over the channels
    for ii = 1:nr_ch
        %% Segmentize the signals
        avg = mean(signal(:,ii)); 
        s = std(signal(:,ii));
        if nr_segm_new > 1
            segmented_ecg = (double(reshape(double(signal(:,ii)),segm_length,nr_segm_new))); 
%             segmented_ecg = segmented_ecg-(repmat(avg,size(segmented_ecg,1),size(segmented_ecg,2)));
%             segmented_ecg = segmented_ecg./(repmat(s,size(segmented_ecg,1),size(segmented_ecg,2)));
        else
            segmented_ecg = (signal(:,ii)-avg)/s;
        end
            
        % Go through each segment
        for iii = 1:nr_segm_new
            %% Make the segments five seconds longer on each side
            if nr_segm_new > 1
                if iii == 1
                    % End the segment five seconds later
                    segm = [segmented_ecg(:,1); segmented_ecg(1:fs*5,2)]; 
                elseif iii==nr_segm_new
                    % Start the segment five seconds earlier
                    segm = [segmented_ecg(end-(fs*5)+1:end,iii-1); segmented_ecg(1:end-(new_signal_length-original_signal_length),iii)]; 
                else
                    % Start the segement five seconds earlier and end five seconds later
                    segm = [segmented_ecg(end-(fs*5)+1:end,iii-1);segmented_ecg(:,iii);segmented_ecg(1:fs*5,iii+1)]; 
                end
            else 
                % Take the whole segment
                segm = segmented_ecg;
            end

            %% Envelope the signal
            % Get the time
            time = (1:length(segm))/fs;

            % Upper envelope
            up_env = env_secant(time,segm,env,'top');

            % Lower envelope
            low_env = env_secant(time,segm,env,'bottom');

            % Enveloped segment
            env_segm = up_env-low_env;

            %% Detect peaks in the enveloped segment
            % Get the frequency ratio with 250 Hz
            freq_ratio = fs/250;
            
            % Get the envelope ratio
            env_ratio = env/(0.15*fs);

            % Peaks are detected by checking whether the derivative is positive and 
            % if the peaks last long enough (step size). The smaller the step 
            % size, the more sure that all R-peaks are detected, but also the 
            % more non-R-peaks are detected. To avoid local minima, a step 
            % (~=1) is used. The original step size is 20 samples = 80ms
            % for an envelope size of 150ms and a sampling frequency of
            % 250 Hz.

            all_peaks = calc_pos_opt(env_segm,floor(20*freq_ratio*env_ratio),1); 
            
%             figure; plot(env_segm); hold on; scatter(all_peaks,env_segm(all_peaks));
 
            % Only continue if enough peaks have been detected
            if length(all_peaks) > time(end)/2
                % Adaptive thresholding
                [R_peak_segm, RR_int_segm, ~] = adaptive_thresholding(all_peaks, env_segm, fs, avgHR, time);
                
%                 figure; plot(env_segm); hold on; scatter(R_peak_segm,env_segm(R_peak_segm));
                
                if postproc && length(RR_int_segm) > 5
                    % Post-process the RR-intervals
                    [R_peak_segm] = post_processing(all_peaks, R_peak_segm, RR_int_segm, env_segm, fs);
                end
            else
                % Set the R and RR variables empty
                R_peak_segm = [];
            end
            
            % Remove excessive R-peaks on both sides
            if length(R_peak_segm) > time(end)/2
                if nr_segm_new > 1
                    if iii == 1
                        % Remove the last five seconds
                        R_peak_segm = R_peak_segm(R_peak_segm<=segm_length);
                    elseif iii == nr_segm_new
                        % Remove the first five seconds
                        R_peak_segm = R_peak_segm(R_peak_segm>fs*5)-(fs*5);
                    else
                        % Remove the first and last five seconds
                        R_peak_segm = R_peak_segm(R_peak_segm>fs*5)-fs*5;
                        R_peak_segm = R_peak_segm(R_peak_segm<=segm_length);
                    end
                end
                R_peak_segm = reshape(R_peak_segm,1,length(R_peak_segm));
            else
                R_peak_segm = [];
            end

             % Store the R-peaks of this channel in a cell
             R_peak_ch{iii} = unique(R_peak_segm);
        end
        
        %% Merge segments
        R_peak_temp = [];
        
        % Loop over the segments
        for iii = 1:nr_segm_new
            R_peak_temp = [R_peak_temp,R_peak_ch{iii}+segm_length*(iii-1)]; %#ok
        end
        
        %% R-peak correction in the ECG signal
        % Define window width (originally this was equal to the envelope size)
        window = round(0.065*fs);
        
        R_peak_temp = peaks_in_ecg(signal(:,ii),R_peak_temp,window);
        
        % Take only the unique R-peaks
        R_peak_temp = unique(R_peak_temp);

        % Get the RR-intervals in ms
        RR_int_temp = 1000*(diff(R_peak_temp)/fs);
            
        if ~isempty(R_peak_temp)
            %% Correct the R-peaks for ectopics
            if ect
                [R_peak_temp,~,~] = ectopic_detection_correction(fs,RR_int_temp,R_peak_temp);
            end
        end
        
        %% Store the R-peaks and RR-intervals
        if ~isempty(R_peak_temp)
            R_peak{ii} = round(unique(R_peak_temp));
            RR_int{ii} = 1000*(diff(R_peak{ii})/fs);
        else
            R_peak{ii} = [];
            RR_int{ii} = [];
        end        
    end
else
    % Set check
    check=0;
end
end




%%%%%% Functions
function [R_peak_segm, RR_int_segm, RR_avg] = adaptive_thresholding(peaks,signal,fs,avgHR,time)

    % Get the standard deviation of the enveloped signal
    STD = std(signal);
    
    % Set the average RR-interval
    RR_avg = zeros(1,length(peaks)-2);
    RR_avg(1) = 60000/avgHR;

%     % Set the average signal and noise peak
%     peak_avg = 0.7*trimmean(signal(peaks),20);
%     noise_avg = peak_avg/2;
    
    peak_avg = 0.7*(median(signal(peaks))+mean(signal(peaks)))/2;
    noise_avg = 0.3*(median(signal(peaks))+mean(signal(peaks)))/2;

    % Set thresholds
    peak_threshold = noise_avg + 0.25*(peak_avg-noise_avg);
    noise_threshold = peak_threshold/2;

    % Pre-allocate
    R_peak_segm = zeros(1,length(peaks));
    RR_int_segm = zeros(1,length(peaks)-1);

    % Select the first R-peak
    count = 1;
    while R_peak_segm(1) == 0
        if signal(peaks(count)) > noise_threshold
            R_peak_segm(1) = peaks(count);
        else
            count = count + 1;
        end
    end

    % Start the peak counter
    peak_counter = 1;

    % Loop over all peaks, starting from the second
    for idx = count + 1 : length(peaks) 

        % Select the peak
        peak = signal(peaks(idx));

        % Define the lower limit of the RR-interval
        RR_lower_limit = 0.65*RR_avg(idx-1);
        if RR_lower_limit < 200
            RR_lower_limit = 200; % Lower limit for RRI is 200 ms --> Includes every condition
        end

        % Define the upper limit of the RR-interval
        RR_upper_limit = 1.35*RR_avg(idx-1);
        if RR_upper_limit > 1600
            RR_upper_limit = 1600; % Upper limit for RRI is 1600 ms
        end

        % Determine whether the peak is a signal or a noise peak
        if peak >= peak_threshold % Peak is an R-peak

            % Make the peak amplitude smaller if it is too high
            % NOTE: This was previously done before the peak
            % search, but this might cause the location of the
            % peak to deviate from the actual location.
            if peak > 4*STD
                peak = 4*STD;
            end

            % Adjust the peak average
            peak_avg = 0.125*peak + 0.875*peak_avg;

            % Add signal peak
            peak_counter = peak_counter + 1;
            R_peak_segm(peak_counter) = peaks(idx);

            % Adjust RR interval
            RR_int_segm(peak_counter-1) = 1000*(R_peak_segm(peak_counter)-R_peak_segm(peak_counter-1))/fs; % In ms

        elseif peak <= noise_threshold % Peak is a noise peak
            % Adjust the noise average
            noise_avg = 0.125*peak + 0.875*noise_avg;

        else % searchback procedure
            if peak_counter > 1

                % Check if the candidate RR-interval is within
                % the boundaries of the lower and upper
                % RR-limits with the previous R-peak
                if time(peaks(idx)) >= RR_lower_limit/1000+time(R_peak_segm(peak_counter))...
                        && time(peaks(idx)) <= RR_upper_limit/1000+time(R_peak_segm(peak_counter))
                    
                    % Check if we are at the end, if the next peak is a
                    % signal peak, if it is a signal peak, check the time
                    % interval
                    if idx == length(peaks) || signal(peaks(idx+1)) < peak_threshold || (signal(peaks(idx+1)) >= peak_threshold && time(peaks(idx+1))-time(peaks(idx)) > 0.2)
                        % Adjust the peak average
                        peak_avg = 0.25*peak + 0.75*peak_avg;

                        % Add R-peak
                        peak_counter = peak_counter + 1;
                        R_peak_segm(peak_counter) = peaks(idx);

                        % Adjust RR interval
                        RR_int_segm(peak_counter-1) = 1000*(R_peak_segm(peak_counter)-R_peak_segm(peak_counter-1))/fs; % In ms
                      
                    else % Peak is a noise peak
                        % Adjust the noise average
                        noise_avg = 0.125*peak+0.875*noise_avg;
                    end

                else % Peak is a noise peak
                    % Adjust the noise average
                    noise_avg = 0.125*peak+0.875*noise_avg;
                end
            else % Peak is a noise peak
                % Adjust the noise average
                noise_avg = 0.125*peak+0.875*noise_avg;
            end
        end

        % Adjust the RR-average
        if peak_counter > 9
            RR_avg(idx) = mean(RR_int_segm(peak_counter-8:peak_counter-1));
        elseif peak_counter > 1
            RR_avg(idx) = 0.5*RR_avg(idx-1)+0.5*median(RR_int_segm(1:peak_counter-1));
        end

        % Adjust thresholds
        peak_threshold = noise_avg + 0.25*(peak_avg-noise_avg);
        noise_threshold = peak_threshold/2;
    end

    % Get the actual R-peaks
    R_peak_segm = R_peak_segm(1:peak_counter);

    % Get the actual RR-intervals, in ms
    RR_int_segm = RR_int_segm(1:peak_counter-1);
    
%     figure; 
%     plot(signal)
%     hold on;
%     scatter(R_peak_segm,signal(R_peak_segm))
end

function R_peak_segm = post_processing(all_peaks, R_peak_segm, RR_int_segm, signal, fs)
    % To check whether RRI is correct. This will be done by comparing
    % the RR-interval with the mean of the 3 past RR-intervals. A past 
    % RR-interval should not be a 'problem' interval.    
    
    % Set basic variables
    a = 0.3; % 30% difference
    d = 0.5; % 50% difference
    max_k = 6; % Maximum ammount of previous RR-intervals

    % Get the amount of RR-intervals
    len = length(RR_int_segm);

    % Pre-allocate
    large_prob = zeros(1,len); % Index of large problems
    RRref_all = zeros(1,len); % Reference RR intervals

    % Set first RR-reference
    RR_ref = trimmean(RR_int_segm(1:5),50);
    RRref_all(1) = RR_ref;

    % Set loop start
    counter = 1;

    % Loop over the different RR-intervals
    while counter < len-1
        % Define the RR-ref if you are further than the first RR-interval
        if counter > 2

            % Set variables
            k = 1; % Number of previous RR-intervals
            sum_ref = 0; % Sum of the good previous RRI's
            num = 0; % Number of good previous RRI's
            wght = [0.5 0.3 0.2]; % Weights: the farther away, the lower the weight

            while counter-k > 0 && num < 3 && k < max_k

                % Check if the previous RRI is a big problem
                if large_prob(counter-k) == 1
                    prob_check = 1;
                elseif large_prob(counter-k) == 0
                    prob_check = 0;
                end

                % If this is not the case, use it for reference
                if ~prob_check
                    num = num + 1;
                    sum_ref = sum_ref + wght(num)*RR_int_segm(counter-k);
                end
                k = k+1;
            end

            % If no good previous RR-intervals are found, this probably means that 
            % this is the new normal, so just take the three previous
            if num <= 1
                k = 1;
                while counter-k > 0 && num < 3 
                    num = num + 1;
                    sum_ref = sum_ref + wght(num)*RR_int_segm(counter-k);
                    k = k+1;
                end
            end

            % Define the reference variables
            RR_ref = sum_ref/sum(wght(1:num));
            no_RR_ref = nnz(RRref_all);
            RRref_all(no_RR_ref+1) = RR_ref;
        end

        % Check for too small RR-intervals
        if RR_int_segm(counter) < (1-a)*RR_ref || RR_int_segm(counter) < 200
            % Set looping variables
            counter1 = 1; 
            counter2 = 0;
            test = 1;

            % Get the sum with the previous RR-interval
            try
                RRnew1 = RR_int_segm(counter) + RR_int_segm(counter-1);
            catch
                RRnew1 = [];
            end

            % Get the sum with the next RR-interval
            sumRR = RR_int_segm(counter) + RR_int_segm(counter + counter1);
            RRnew2 = sumRR;

            % Loop until a good summation is found
            while test == 1 

                % Only proceed if we are not at the end of the signal 
                if counter + counter1 < len
                    % If the new RR-interval is still too small
                    if sumRR < (1-a)*RR_ref || sumRR < 200 
                        % Adjust counter
                        counter1 = counter1 + 1;

                        % Add the next RR-interval
                        sumRR = sumRR + RR_int_segm(counter+counter1);

                    % Good summation   
                    elseif abs(sumRR-RR_ref) < a*RR_ref
                        % Set new looping variable                     
                        test2 = 1;

                        % Adjust the second counter
                        counter2 = counter2 + 1;

                        % Loop until the best summation is found
                        while test2 == 1

                            % Only proceed if we are not at the end of the
                            % signal
                            if counter + counter1 + counter2 < len

                                % Add the next RR-interval
                                sumRR2 = sumRR + RR_int_segm(counter + counter1 + counter2);

                                % Compute the difference of both sums with the
                                % RR-reference
                                diff1 = abs(RR_ref-sumRR) / RR_ref;
                                diff2 = abs(RR_ref-sumRR2) / RR_ref;

                                % If the initial sum is best
                                if diff1 < diff2 
                                    % Stop both loops
                                    test = 0; 
                                    test2 = 0; 

                                    % Store the initial sum as new RR-interval
                                    RRnew2 = sumRR;

                                    % Adjust the second counter (necessary for
                                    % later)
                                    counter2 = counter2-1;

                                % If the second sum is best    
                                else 
                                    % Redefine the initial sum and go through
                                    % the loop again
                                    sumRR = sumRR2;

                                    % Adjust the counter
                                    counter2 = counter2+1;

                                    % Stop the loop if we are at the end of the
                                    % signal
                                    if counter + counter1 + counter2 > len
                                        test = 0;
                                        test2 = 0;
                                        RRnew2 = sumRR2; 
                                    end
                                end
                            else
                                % Stop the loop
                                test = 0;
                                test2 = 0;
                                RRnew2 = sumRR;
                            end
                        end
                    % Too big    
                    else 
                        % Stop the loop
                        test = 0;

                        % Select the previous sum as correct
                        RRnew2 = sumRR - RR_int_segm(counter + counter1);

                        % Adjust the first counter (necessary for later)
                        counter1 = counter1 - 1;
                    end
                else
                    % Stop loop
                    test = 0;
                    RRnew2 = sumRR;
                end
            end

            % If no better option has been found,
            % or if the new option is better than the
            % combination of the previous RR-intervals

            if ~isempty(RRnew1) && abs(RRnew1 - RR_ref) < abs(RRnew2 - RR_ref) && abs(RRnew1-RR_ref) < a*RR_ref
                % Set RRnew
                RRnew = RRnew1;

                % Replace RR-interval
                RR_int_segm(counter - 1) = RRnew1;

                % Change RR
                RR_int_segm(counter) = [];

                % Change large prob vector (has equal length as RR)
                large_prob(counter) = [];

                % Change RRref_all (is one smaller than RR)
                if counter + 1 < len-1
                    RRref_all(counter) = [];
                end

                % Get new length
                len = length(RR_int_segm); 
                
                % Make sure that the new interval is checked
                counter = counter-1;
            else
                % Set RRnew
                RRnew = RRnew2;

                % Replace RR-interval
                RR_int_segm(counter) = RRnew;

                % Change the variables depending on circumstances
                if counter1 >= 1
                    % Change RR
                    RR_int_segm(counter + 1 : counter + counter1 + counter2) = [];

                    % Change large prob vector (has equal length as RR)
                    large_prob(counter + 1 : counter + counter1 + counter2) = [];

                    % Change RRref_all (is one smaller than RR)
                    if counter+1 < len-1
                        RRref_all(counter+1:min([counter+counter1+counter2, len-1])) = [];
                    end

                    % Get new length
                    len = length(RR_int_segm); 
                end
            end


            % Check if the interval is a big problem
            if abs(RRnew-RR_ref)/RR_ref > d || RRnew < 200 % Check for big problems
                large_prob(counter) = 1;
            end

        % Check for too big RR-intervals   
        elseif RR_int_segm(counter) > (1+a)*RR_ref || RR_int_segm(counter) > 1600

            % Convert RR-intervals to samples
            RR_samples = RR_int_segm*fs/1000;

            % Define R-peak positions
            R_peak_positions = cumsum(RR_samples(1:counter)) + R_peak_segm(1);

            % Define beginning and end of the interval
            if counter > 1
                begin_int = R_peak_positions(end-1);
            else
                begin_int = R_peak_segm(1);
            end
            end_int = R_peak_positions(end);

            % Get the first peak that is bigger than the beginning of the
            % interval
            test_peak = all_peaks(all_peaks > begin_int);
            test_peak = test_peak(test_peak < end_int);

            if length(test_peak) > 1
                % Get the test RR interval
                RR_test = 1000*(test_peak-begin_int)/fs;

                % Find the smallest difference with the
                % reference RR interval
                [~,I] = min(abs(RR_test-RR_ref));
                test_peak = test_peak(I);
            end

            % See if that peak is before the end of the
            % interval and not too small
            if ~isempty(test_peak) &&...
                    signal(test_peak)                 > prctile(signal,50)          && ...
                    (((test_peak-begin_int)/fs)*1000) > (1-a)*RR_ref                &&...
                    (((test_peak-begin_int)/fs)*1000) > 200                         &&...
                    (((end_int-test_peak)/fs)*1000)   > (1-a)*RR_ref                && ...
                    (((end_int-test_peak)/fs)*1000)   > 200

                % Add an RR interval
                if counter > 1
                    RR_int_segm = [RR_int_segm(1:counter-1) (((test_peak-begin_int)/fs)*1000) (((end_int-test_peak)/fs)*1000) RR_int_segm(counter+1:end)];
                elseif counter+1 > len
                    RR_int_segm = [RR_int_segm(1:counter-1) (((test_peak-begin_int)/fs)*1000) (((end_int-test_peak)/fs)*1000)];
                else
                    RR_int_segm = [(((test_peak-begin_int)/fs)*1000) (((end_int-test_peak)/fs)*1000) RR_int_segm(counter+1:end)];
                end

                % Get new length
                len = length(RR_int_segm);

                % Adjust the large problem variable
                large_prob = [large_prob 0]; %#ok

                % Adjust the RR-ref all variable
                RRref_all = [RRref_all 0]; %#ok

                % Make sure that the new interval is checked
                counter = counter-1;

            else
                % State that the interval is a large problem and should not
                % be included for the reference intervals
                large_prob(counter) = 1;
            end
        end

        % Go one peak further
        counter = counter + 1;
    end

    % Define the R-peak positions 
    RR_int_segm = RR_int_segm*fs/1000; % Go back to samples
    R_peak_segm(2:length(RR_int_segm)+1) = round(cumsum(RR_int_segm) + R_peak_segm(1));  
end  


function R_peak_temp = peaks_in_ecg(signal,R_peak_temp,window)
    
%     % Get the percentage of positive R-peaks
%     perc_pos = sum(sign(signal(R_peak_temp)))/length(R_peak_temp);

    for idx = 1:length(R_peak_temp) 
        % Get the samples from the R-peak minus the envelope length
        % to the R-peak

        p = signal(max(1, R_peak_temp(idx) - window) : R_peak_temp(idx)); 

        % Do the same but for the time locations
        t = max(1,R_peak_temp(idx)-window) : R_peak_temp(idx);

        % Search for the maximum in p
        [~,ip] = max(p);
        
%         % Search for the extremum in p
%         if perc_pos >= 0.5
%             [~,ip] = max(p); 
%         else
%             [~,ip] = min(p);
%         end

        % Adjust the R-peak
        try
            if signal(t(ip(end))) >= signal(t(ip(end))-1) && signal(t(ip(end))) >= signal(t(ip(end))+1)
                R_peak_temp(idx) = t(ip(end)); 
            end
        catch
            R_peak_temp(idx) = t(ip(end));
        end
    end
end
