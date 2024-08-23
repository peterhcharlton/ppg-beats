function [peaks, onsets] = mmpdv2_beat_detector(sig,fs)
% MMPDV2_BEAT_DETECTOR  PDA PPG beat detector.
%   MMPDV2_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
%   using the 'Mountaineer's Method for Peak Detection' (v.2) beat detector
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
%   E. J. A. Prada, 'The mountaineer's method for peak detection in photoplethysmographic signals,' Revista Facultad de Ingenieria, vol. 90, pp. 42-50, 2019. <https://doi.org/10.17533/udea.redin.n90a06>
%   
%   # Author
%   * E. J. Arguello Prada - wrote the code: 'Alpinista_simple_4_todos'
%   * Peter H. Charlton - did very little, just wrote this wrapper
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Version
%   1.0
%   
%   # License - MIT
%      Copyright (c) 2024 E. J. Arguello Prada and Peter H. Charlton
%      Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
%      The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
%      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

[~, peaks] = Alpinista_simple_4_todos(sig,fs);

peaks = tidy_beats(peaks);

onsets = pulse_onsets_from_peaks(sig, peaks);

end
% The mountaineer's method for PPG peaks detection
% a = PPI = Peak-to-peak intervals
% b = PWV = Pulse width or valley-to-valley intervals
% c = Tc = Crest time
% d = Td = Downslope time
% e = PPGA = Systolic peak amplitude
% f = PPGA_abs = Pulse amplitude
% g = PPGA_norm = Normalized pulse amplitude
% h = PPAT = Pulse area
% k = PPA = Amplitude differences between peaks
% l = PPB = Diferencia entre mínimos

function [PPI, t_pico, PWV, Tc, Td, PPGA, PPGA_abs, PPGA_norm, PPAT, PPA, PPB] = Alpinista_simple_4_todos(senal,fsample)   % adjusted by PC to return "t_pico" (beat timings) too
% Initialization 
cuenta_ascenso = 0;
ascenso_max = round(0.6*fsample*0.15);            % Threshold depends on the sampling rate.
primer_pulso = false;

PPI_value = 1.0;
refractory = 0.35;       % For non-pediatric PPG recordings, set this value to 0.65
j=1;

% Auxiliar variables
contador(1) = 0;
umbral(1) = 0;
intervalo(1) = PPI_value;
refract(1) = refractory;


y1 = senal;
num_muestras = length(y1);
tiempo = (1:num_muestras);
% Execution time measurement using tic-toc function
% tic  % commented by PC
for i=2:num_muestras       
    % Detection starts here!
    if y1(i) > y1(i-1)
        cuenta_ascenso = cuenta_ascenso + 1;
    else
        if cuenta_ascenso >= ascenso_max
            % If the first peak has not been detected...
            if primer_pulso == false
                primer_pulso = true;
                amp_pico(j) = y1(i-1);
                t_pico(j) = tiempo(i-1);
                
                % Valley detection
                amp_min(j) = y1(i - (1 + cuenta_ascenso));
                t_min(j) = tiempo(i - (1 + cuenta_ascenso));
                
                ascenso_max = round(0.6*fsample*0.15); 
                j = j + 1;
            else
                % If the first peak has been detected, the next one must be
                % far enough to the first one
                if (tiempo(i)-t_pico(j-1))/fsample > 1.2*PPI_value || cuenta_ascenso > round((1.75*ascenso_max)/(0.6)) % Esta condición hace la diferencia!   
                    amp_pico(j) = y1(i-1);
                    t_pico(j) = tiempo(i-1);
                    
                    % Valley detection
                    amp_min(j) = y1(i - (1 + cuenta_ascenso));
                    t_min(j) = tiempo(i - (1 + cuenta_ascenso));
                    
                    ascenso_max = round(0.6*fsample*0.15);
                    PPI(j-1) = (t_pico(j) - t_pico(j-1))/fsample;
                    
                    % Fiducial-based feature calculation (see above)
                    PPA(j-1) = amp_pico(j) - amp_pico(j-1);                 % Variabilidad sistólica o diferencias entre máximos
                    PPB(j-1) = amp_min(j) - amp_min(j-1);                   % Variabilidad diastólica o diferencias entre mínimos
                    PWV(j-1) = (t_min(j) - t_min(j-1))/fsample;             % Ancho del pulso
                    PPAT(j-1) = ((t_min(j)-t_min(j-1))/fsample)*(amp_pico(j-1)-amp_min(j-1))/2;   % Área del pulso
                    
                    % Peak-to-peak interval and refractory period are
                    % updated
                    PPI_value = PPI(j-1);
                    refractory = 0.35;
%                     primer_pulso = false;
                    num_picos = j;
                    j = j + 1;
                else
                    if (tiempo(i)-t_pico(j-1))/fsample > refractory
                        amp_pico(j) = y1(i-1);
                        t_pico(j) = tiempo(i-1);
                        % Valley detection
                        amp_min(j) = y1(i - (1 + cuenta_ascenso));
                        t_min(j) = tiempo(i - (1 + cuenta_ascenso));
                        
                        ascenso_max = round(0.6*cuenta_ascenso); % Se actualiza el umbral...
                        PPI(j-1) = (t_pico(j) - t_pico(j-1))/fsample;
                        
                        % Fiducial-based feature calculation (see above)
                        PPA(j-1) = amp_pico(j) - amp_pico(j-1);                 % Variabilidad sistólica o diferencias entre máximos
                        PPB(j-1) = amp_min(j) - amp_min(j-1);                   % Variabilidad diastólica o diferencias entre mínimos
                        PWV(j-1) = (t_min(j) - t_min(j-1))/fsample;             % Ancho del pulso
                        PPAT(j-1) = ((t_min(j)-t_min(j-1))/fsample)*(amp_pico(j-1)-amp_min(j-1))/2;   % Área del pulso
                        
                        % Peak-to-peak interval and refractory period are
                        % updated
                        PPI_value = PPI(j-1);
                        refractory = 0.75*PPI_value;
                        num_picos = j;
                        j = j + 1;
                    end                        
                end             
            end
        end        
        cuenta_ascenso = 0;
    end
    contador(i) = cuenta_ascenso;
    umbral(i) = ascenso_max;
    intervalo(i) = PPI_value;
    refract(i) = refractory;
end
%toc  % commented by PC

% One-pulse fiducial-based feature calculation
for j = 1:num_picos
    if length(t_pico) ~= length(t_min)
        if length(t_pico) > length(t_min)
            t_pico = t_pico(1:end-(length(t_pico) - length(t_min)));
        else
            t_min = t_min(1:end-(length(t_min) - length(t_pico)));
        end
    end
    Tc(j) = (t_pico(j)-t_min(j))/fsample;
    PPGA(j) = amp_pico(j);
    PPGA_abs(j) = PPGA(j) - amp_min(j);
    PPGA_norm(j) = PPGA_abs(j)/PPGA(j);
    if j < num_picos
        Td(j) = (t_min(j+1) - t_pico(j))/fsample;
    end
end


% Plot the detection results
do_plot = 0; % added by PC
if do_plot % added by PC
    figure
    subplot(2,1,1)
    plot(tiempo./fsample,y1,t_pico./fsample,amp_pico,'or','MarkerSize',12);
    hold on;
    plot(t_min./fsample,amp_min,'dk','MarkerSize',12);                           % ... y con esto se marcan los mínimos
    % plot(tiempo./fsample,contador,':k',tiempo./fsample,umbral,'r');
    % plot(tiempo./fsample,intervalo,'r',tiempo./fsample,refract,':k');
    hold off;
    xlabel('Time (s)','fontsize',14);
    ylabel('Magnitude (V)','fontsize',14);
    xlim([0 length(y1)/fsample]);
    set(gca,'fontsize',12);     % Permite modificar el tamaño de la fuente para los valores en los ejes
    set(gcf,'color','w');
end % added by PC

end
