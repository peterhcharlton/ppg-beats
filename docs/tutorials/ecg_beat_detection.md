# ECG Beat Detection

This page provides examples of how to use the toolbox for ECG beat detection.

## Detecting beats in an ECG signal

This tutorial demonstrates how to detect beats in an ECG signal using algorithms in the toolbox.

- Download ECG data from ???
- Load this data file into Matlab. The file contains a single variable named _S_, which is a structure containing a 1-minute ECG signal (with the ECG signal given in two fields: _v_ contains the ECG samples, and _fs_ is the sampling frequency in Hz).
- Use the following Matlab commands to detect beats in the ECG, and then plot the ECG signal and detected beats:

```matlab
options.win_durn = 20;     % assess agreement between beat detectors using 20-second windows
[beat_inds, qual] = detect_ecg_beats(S.v, S.fs, options);     % detect beats and assess quality of beat detections

figure('Position', [20,20,1000,350])     % Setup figure
subplot('Position', [0.05,0.17,0.92,0.82])
t = [0:length(S.v)-1]/S.fs;             % Make time vector
plot(t, S.v, 'b'), hold on,             % Plot ECG signal
plot(t(find(~qual)), S.v(~qual), 'k')         % Plot low quality segments
plot(t(beat_inds), S.v(beat_inds), 'or'),       % Plot detected beats
ftsize = 20;                            % Tidy up plot
set(gca, 'FontSize', ftsize, 'YTick', []);
ylabel('ECG', 'FontSize', ftsize),
xlabel('Time (s)', 'FontSize', ftsize)
legend({'High quality', 'Low quality', 'Beats'}, 'FontSize', ftsize, 'Location', 'northoutside','NumColumns',3)     % Make plot legend
```

This results in the following:

... IMAGE AND DESCRIPTION TO BE INSERTED ...


## Assessing the quality of beat detections