function [peaks, onsets] = pwd_beat_detector(sig, fs)
% PWD_BEAT_DETECTOR  PWD PPG beat detector.
%   PWD_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
%   using the 'Pulse Wave Delineator' beat detector
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
%   B. N. Li et al., 'On an automatic delineator for arterial blood pressure waveforms,' Biomedical Signal Processing and Control, vol. 5, no. 1, pp. 76-81, 2010. <https://doi.org/10.1016/j.bspc.2009.06.002>
%   
%   # Author
%   * LI Bing Nan - created the original code
%   * Peter H. Charlton - did very little, just wrote this wrapper
%   
%   # Documentation
%   <https://ppg-beats.readthedocs.io/>
%   
%   # Version
%   1.0
%   
%   # Source
%   <https://www.mathworks.com/matlabcentral/fileexchange/29484-pulse-waveform-delineator>
%
%   # Licence
%      The 2-Clause BSD License (see the copyright notice, list of conditions, and disclaimer below)

[onsets,peaks,~] = delineator(sig,fs);
onsets = onsets(:);
peaks = peaks(:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This section is 'Pulse Waveform Delineator' code
%
% from: https://www.mathworks.com/matlabcentral/fileexchange/29484-pulse-waveform-delineator
%
% Copyright (c) 2010, LI Bing Nan
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [onsetp,peakp,dicron] = delineator(abpsig,abpfreq)
% This program is intended to delineate the fiducial points of pulse waveforms
% Inputs:
%   abpsig: input as original pulse wave signals;
%   abpfreq: input as the sampling frequency;
% Outputs:
%   onsetp: output fiducial points as the beginning of each beat;
%   peakp: output fiducial points as systolic peaks;
%   dicron: output fiducial points as dicrotic notches;
% Its delineation is based on the self-adaptation in pulse waveforms, but
% not in the differentials.
% Reference:
%   BN Li, MC Dong & MI Vai (2010) 
%   On an automatic delineator for arterial blood pressure waveforms
%   Biomedical Signal Processing and Control 5(1) 76-81.
% LI Bing Nan @ University of Macau, Feb 2007
%   Revision 2.0.5, Apr 2009
%Initialization
peakIndex=0;
onsetIndex=0;
dicroIndex=0;
stepWin=2*abpfreq;
closeWin=floor(0.1*abpfreq);    %invalide for pulse beat > 200BPM
sigLen=length(abpsig);
peakp=[];
onsetp=[];
dicron=[];
%lowpass filter at first
coh=25;                     %cutoff frequency is 25Hz
coh=coh*2/abpfreq;
od=3;                       %3rd order bessel filter
[B,A]=besself(od,coh);
abpsig=filter(B,A,abpsig);
abpsig=10*abpsig;
abpsig=smooth(abpsig);
%Compute differentials
ttp=diff(abpsig);
diff1(2:sigLen)=ttp;
diff1(1)=diff1(2);
diff1=100*diff1;
clear ttp;
diff1=smooth(diff1);
if sigLen>12*abpfreq
    tk=10;
elseif sigLen>7*abpfreq
    tk=5;
elseif sigLen>4*abpfreq
    tk=2;
else
    tk=1;
end
%Seek avaerage threshold in original signal
if tk>1             %self-learning threshold with interval sampling
    tatom=floor(sigLen/(tk+2));
    for ji=1:tk       %search the slopes of abp waveforms
        sigIndex=ji*tatom;
        tempIndex=sigIndex+abpfreq;
        [tempMin,jk,tempMax,jl]=seeklocales(abpsig,sigIndex,tempIndex);
        tempTH(ji)=tempMax-tempMin;
    end
    abpMaxTH=mean(tempTH);
else
    [tempMin,jk,tempMax,jl]=seeklocales(abpsig,closeWin,sigLen);
    abpMaxTH=tempMax-tempMin;
end
clear j*;
clear t*;
abpMaxLT=0.4*abpMaxTH;
%Seek pulse beats by MinMax method
% diffIndex=1;
diffIndex=closeWin;             %Avoid filter distortion
while diffIndex<sigLen
    tempMin=abpsig(diffIndex);   %Initialization
    tempMax=abpsig(diffIndex);
    tempIndex=diffIndex;
    tpeakp=diffIndex;        %Avoid initial error
    tonsetp=diffIndex;      %Avoid initial error
    while tempIndex<sigLen
        %If no pulses within 2s, then adjust threshold and retry
        if (tempIndex-diffIndex)>stepWin
%             tempIndex=diffIndex-closeWin;
            tempIndex=diffIndex;
            abpMaxTH=0.6*abpMaxTH;
            if abpMaxTH<=abpMaxLT
                abpMaxTH=2.5*abpMaxLT;
            end
            break;
        end
        if (diff1(tempIndex-1)*diff1(tempIndex+1))<=0  %Candidate fiducial points
            if (tempIndex+5)<=sigLen
                jk=tempIndex+5;
            else
                jk=sigLen;
            end
            if (tempIndex-5)>=1
                jj=tempIndex-5;
            else
                jj=1;
            end
            %Artifacts of oversaturated or signal loss?
            if (jk-tempIndex)>=5
                for ttk=tempIndex:jk
                    if diff1(ttk)~=0
                        break;
                    end
                end
                if ttk==jk
                    break;          %Confirm artifacts
                end
            end
            if diff1(jj)<0          %Candidate onset
                if diff1(jk)>0
                    [tempMini,tmin,ta,tb]=seeklocales(abpsig,jj,jk);
                    if abs(tmin-tempIndex)<=2
                        tempMin=tempMini;
                        tonsetp=tmin;
                    end
                end
            elseif diff1(jj)>0      %Candidate peak
                if diff1(jk)<0
                    [tc,td,tempMaxi,tmax]=seeklocales(abpsig,jj,jk);
                    if abs(tmax-tempIndex)<=2
                        tempMax=tempMaxi;
                        tpeakp=tmax;
                    end
                end
            end
            if ((tempMax-tempMin)>0.4*abpMaxTH)   %evaluation
                if ((tempMax-tempMin)<2*abpMaxTH)
                    if tpeakp>tonsetp
                        %If more zero-crossing points, further refine!
                        ttempMin=abpsig(tonsetp);
                        ttonsetp=tonsetp;
                        for ttk=tpeakp:-1:(tonsetp+1)
                            if abpsig(ttk)<ttempMin
                                ttempMin=abpsig(ttk);
                                ttonsetp=ttk;
                            end
                        end
                        tempMin=ttempMin;
                        tonsetp=ttonsetp;
                            
                        if peakIndex>0
                            %If pulse period less than eyeclose, then artifact
                            if (tonsetp-peakp(peakIndex))<(3*closeWin)
                                %too many fiducial points, then reset
                                tempIndex=diffIndex;                                
                                abpMaxTH=2.5*abpMaxLT;
                                break;
                            end
                            
                            %If pulse period bigger than 2s, then artifact
                            if (tpeakp-peakp(peakIndex))>stepWin
                                peakIndex=peakIndex-1;
                                onsetIndex=onsetIndex-1;
                                if dicroIndex>0
                                    dicroIndex=dicroIndex-1;
                                end
                            end
                            if peakIndex>0
                                %new pulse beat
                                peakIndex=peakIndex+1;
                                peakp(peakIndex)=tpeakp;
                                onsetIndex=onsetIndex+1;
                                onsetp(onsetIndex)=tonsetp;
                                tf=onsetp(peakIndex)-onsetp(peakIndex-1);
                                to=floor(abpfreq./20);   %50ms
                                tff=floor(0.1*tf);
                                if tff<to
                                    to=tff;
                                end
                                to=peakp(peakIndex-1)+to;
                                te=floor(abpfreq./2);   %500ms
                                tff=floor(0.5*tf);
                                if tff<te
                                    te=tff;
                                end
                                te=peakp(peakIndex-1)+te;
                                tff=seekdicrotic(diff1(to:te));
                                if tff==0
                                    tff=te-peakp(peakIndex-1);
                                    tff=floor(tff/3);
                                end
                                dicroIndex=dicroIndex+1;
                                dicron(dicroIndex)=to+tff;
                                tempIndex=tempIndex+closeWin;
                                break;
                            end
                        end
                        
                        if  peakIndex==0   %new pulse beat
                            peakIndex=peakIndex+1;
                            peakp(peakIndex)=tpeakp;
                            onsetIndex=onsetIndex+1;
                            onsetp(onsetIndex)=tonsetp;
                            tempIndex=tempIndex+closeWin;
                            break;
                        end
                    end
                end
            end
        end
        tempIndex=tempIndex+1;      %step forward
    end
%     diffIndex=tempIndex+closeWin;    %for a new beat
    diffIndex=tempIndex+1;
end
if isempty(peakp),return;end
%Compensate the offsets of lowpass filter
sigLen=length(peakp);
for diffIndex=1:sigLen          %avoid edge effect
    tempp(diffIndex)=peakp(diffIndex)-od;
end
ttk=tempp(1);
if ttk<=0
    tempp(1)=1;
end 
clear peakp;
peakp=tempp;
clear tempp;
sigLen=length(onsetp);
for diffIndex=1:sigLen
    tempp(diffIndex)=onsetp(diffIndex)-od;
end
ttk=tempp(1);
if ttk<=0
    tempp(1)=1;
end 
clear onsetp;
onsetp=tempp;
clear tempp;
if isempty(dicron),return;end
sigLen=length(dicron);
for diffIndex=1:sigLen
    if dicron(diffIndex)~=0
        tempp(diffIndex)=dicron(diffIndex)-od;
    else
        tempp(diffIndex)=0;
    end
end
clear dicron;
dicron=tempp;
clear tempp;

end

function [mini,minip,maxi,maxip]=seeklocales(tempsig,tempbegin,tempend)
tempMin=tempsig(tempbegin);
tempMax=tempsig(tempbegin);
minip=tempbegin;
maxip=tempbegin;
for j=tempbegin:tempend
    if tempsig(j)>tempMax
        tempMax=tempsig(j);
        maxip=j;
    elseif tempsig(j)<tempMin
        tempMin=tempsig(j);
        minip=j;
    end
end
mini=tempMin;
maxi=tempMax;

end

function [dicron]=seekdicrotic(tempdiff)
izcMin=0;
izcMax=0;
itemp=3;
tempLen=length(tempdiff)-3;
dicron=0;
tempdiff=smooth(tempdiff);
while itemp<=tempLen
    if (tempdiff(itemp)*tempdiff(itemp+1))<=0
        if tempdiff(itemp-2)<0
            if tempdiff(itemp+2)>=0
                izcMin=izcMin+1;
                tzcMin(izcMin)=itemp;
            end
        end
%         if tempdiff(itemp-2)>0
%             if tempdiff(itemp+2)<=0
%                 izcMax=izcMax+1;
%                 tzcMax(izcMax)=itemp;
%             end
%         end
    end
    itemp=itemp+1;
end
if izcMin==0     %big inflection
    itemp=3;
    tempMin=tempdiff(itemp);
    itempMin=itemp;
    
    while itemp<tempLen
        if tempdiff(itemp)<tempMin
            tempMin=tempdiff(itemp);
            itempMin=itemp;
        end
        itemp=itemp+1;
    end
    itemp=itempMin+1;
    while itemp<tempLen
        if tempdiff(itemp+1)<=tempdiff(itemp-1)
            dicron=itemp;
            return;
        end
        itemp=itemp+1;
    end
elseif izcMin==1
    dicron=tzcMin(izcMin);
    return;
else
    itemp=tzcMin(1);
    tempMax=tempdiff(itemp);
    itempMax=itemp;
    
    while itemp<tempLen
        if tempdiff(itemp)>tempMax
            tempMax=tempdiff(itemp);
            itempMax=itemp;
        end
        itemp=itemp+1;
    end
    for itemp=izcMin:-1:1
        if tzcMin(itemp)<itempMax
            dicron=tzcMin(itemp);
            return;
        end
    end
end

end

function [diap]=seekdiap(tempabp)
diap=0;
[tt,ti]=max(tempabp);
if ti==0
    diap=floor(length(tempabp)./2);
else
    diap=ti;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%