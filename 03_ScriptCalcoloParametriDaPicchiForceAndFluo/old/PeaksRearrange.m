%rearrange start and end itme points
%Sig = filtered signal
function [StartStopPeak_assessed ModArray] = PeaksRearrange(Sig,Inf,Fs)

%num peaks
NumPeak = size(Inf,1);
%peak duration [points]
PeakDur = round(abs(Inf(:,2)*Fs));
%start e stop peak index
StartStopPeak = [round(abs(Inf(:,1)*Fs))  round(abs(Inf(:,1)*Fs))+PeakDur-1];


StartStopPeak_assessed = StartStopPeak*0;
ModArray = zeros(length(Sig),1);

%
figure
plot(Sig)
%


%trova inizio picco e fine picco (controllando i minimi prima e dopo il picco trovato)
for i=1:NumPeak
   
    %
    hold on
    plot(StartStopPeak(i,1):StartStopPeak(i,2),Sig(StartStopPeak(i,1):StartStopPeak(i,2)),'r')
    %
    
    %counter
    d = i/NumPeak * 100
    
    %first peak
    if i==1
        
        %parte segnale di interesse dove trovare primo minimo subito precedente al primo picco
        indexPrePeak  = 1:StartStopPeak(i,1);
        
        %if only a peak
        if NumPeak==1
            indexPostPeak = StartStopPeak(i,2):length(Sig);
        else
            indexPostPeak = StartStopPeak(i,2):StartStopPeak(i+1,1);
        end
        
    %last peak
    elseif i==NumPeak
        
        %parte segnale di interesse dove trovare primo minimo subito precedente all'ultimo picco e subito dopo
        indexPrePeak   = StartStopPeak(i-1,2):StartStopPeak(i,1);
        indexPostPeak  = StartStopPeak(i,2):length(Sig);
        
    %otherwise
    else
        
        indexPrePeak   = StartStopPeak(i-1,2):StartStopPeak(i,1);
        indexPostPeak  = StartStopPeak(i,2):StartStopPeak(i+1,1);
        
    end
    
    partPrePeak   = Sig(indexPrePeak);
    partPostPeak  = Sig(indexPostPeak);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %peak start
    if length(partPrePeak)>=3
        [Peaks PeaksLocs] = findpeaks( partPrePeak(end)-partPrePeak );
        if isempty(PeaksLocs)
            StartPeak = round(median(indexPrePeak));
        else
            StartPeak = indexPrePeak(PeaksLocs(end));
        end
    else
        StartPeak  = StartStopPeak(i,1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %peak end
    if length(partPostPeak)>=3
        [Peaks PeaksLocs] = findpeaks( partPostPeak(1)- partPostPeak);
        if isempty(PeaksLocs)
            StopPeak = round(median(indexPostPeak));
        else
            StopPeak = indexPostPeak(PeaksLocs(1));
        end
    else
        StopPeak  = StartStopPeak(i,2);
    end
    
    %update
    %%%%%%%%%%%%%%%%%%%%%%%%
    StartStopPeak(i,1) = StartPeak;
    StartStopPeak(i,2) = StopPeak;    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    StartStopPeak_assessed(i,1) = StartPeak;
    StartStopPeak_assessed(i,2) = StopPeak;
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    ModArray(StartPeak:StopPeak) = 1;
    
    
    %
    hold on
    plot(StartStopPeak_assessed(i,1):StartStopPeak_assessed(i,2),Sig(StartStopPeak_assessed(i,1):StartStopPeak_assessed(i,2)),'g')
    %
    
end


end