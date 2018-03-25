%Script To compute Parameters from Force Peaks
%ComputeForceParScript

% close all
clc

%rearrange -> il rearrange non va mai fatto
Rearrange = 0;

if isfield(dataGCamp,'PeaksPar')
    display('Computation of parameters already performed')
else
    
            
    %sampling frequency
    Fs = dataGCamp.Info.Fs;
    %task status
    Status = dataGCamp.status;
    %number of the trials (vector)
    TrialsNum = dataGCamp.InfoTrial.TrialsVector;
    %resting (no force/fluo peaks) interval -> used to quantify f0 for the normalization of fluo
    RestInterval = dataGCamp.Info.IntervalToFluoForceMean;
    
    %num Fluorescence ROIs + Force Signal
%     NumSignals = size(dataGCamp.fluoROI,2)+1;
    NumSignals = 1;
        
    %%%%% Find Peaks  %%%%%%%%%  %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%%
    for iNS=1:NumSignals  
        
        if iNS == 1
            %Force is now positive
            SigOriginal     = -dataGCamp.fx;
            Sig = SigOriginal;
            %Threshold 3std
            Threshold  = mean(Sig(RestInterval(1):RestInterval(2)))+ std(Sig(RestInterval(1):RestInterval(2)))*3;
            
            %Speed
            Speed      = -dataGCamp.speed;
            Threshold_Speed  = mean(Speed)+std(Speed)*3;
            
        else
%             %Fluo ROI
%             SigOriginal     = dataGCamp.fluoROI(:,iNS-1);
%             %filtering     
%             Sig = cheb2LPfilt(SigOriginal,1,2,Fs);
%             SigOriginal = cheb2LPfilt(SigOriginal,6,5,Fs);
%             %Threshold 3std
%             Threshold  = std(Sig(RestInterval(1):RestInterval(2)))*3;
        end
        
        
        %find peaks
        [Inf] = PeaksFinder_v3(Sig ,Threshold, Fs,1);
        
        %if peaks have been found
        if ~isempty(Inf) %if peak
            
            if Rearrange==1
                %rearrange start and end peaks
                [StartStopPeak PeakArray] = PeaksRearrange_v2(Sig,Inf,Fs,SigOriginal,'on');
            else
                PeakDur = round(abs(Inf(:,2)*Fs));
                StartStopPeak = [round(abs(Inf(:,1)*Fs))  round(abs(Inf(:,1)*Fs))+PeakDur-1];
                
                PeakArray = zeros(length(Sig),1);
                for ai=1:size(StartStopPeak,1)
                    PeakArray(StartStopPeak(ai):StartStopPeak(ai)) = 1;
                end
            end
            
            %force rearrangement based on speed peaks
            if iNS == 1         %find peaks
                %speed
                [Inf_Speed] = PeaksFinder_v3(Speed,Threshold_Speed, Fs,1);
                %rearrange start and end peaks of speed
                [StartStopPeak_Speed PeakArray_Speed] = PeaksRearrange_v2(Speed,Inf_Speed,Fs,Speed,'on');
                
                if ~isempty(Inf_Speed) %if peak
                    [StartStopPeak_Force PeakArray_Force] = ForcePeaksCleaner_FromSpeed(StartStopPeak, StartStopPeak_Speed,Fs,PeakArray, Sig, Speed,'off');
                    StartStopPeak = StartStopPeak_Force;
                    PeakArray     = PeakArray_Force;
                else
                    error('no speed peaks')
                end
            end
            
            %%%%% Find Info Peaks and Parameters %%%%%%%%%  %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%%
            
            % parameters from peaks
            [PeaksPar] = ParFromPeaks(SigOriginal, Fs,StartStopPeak, TrialsNum, Status,1);
            
            %store in a new field of dataGCamp
            dataGCamp.PeakArray(:,iNS) = PeakArray;
            dataGCamp.PeaksPar{1,iNS}  = PeaksPar;
            
        end %end if peak
        
    end
end

            
%check force - fluo
if 0
    %ForceCol
    ForceCol = 1;
    %Status to check = [3 4]
    StatusCheck= [3 4]; 
    %to add the struct field: CorrespondanceForceFluoPeaks -> 
    dataGCamp = PeaksFindFluoForceCorresp(dataGCamp,ForceCol, StatusCheck);
end


if 1
    figure
    t     = dataGCamp.t;
    fx    = dataGCamp.PeakArray(:,1);
%     fluo1 = dataGCamp.PeakArray(:,2);
%     fluo2 = dataGCamp.PeakArray(:,3);
%     fluo3 = dataGCamp.PeakArray(:,4);
%     fluo4 = dataGCamp.PeakArray(:,5);
%     fluo5 = dataGCamp.PeakArray(:,6);
%     fluo6 = dataGCamp.PeakArray(:,7);
%     fluo7 = dataGCamp.PeakArray(:,8);
    
    
    plot(t,fx,'r')
    hold on
%     plot(t,fluo1+1.5)
%     plot(t,fluo2+3)
%     plot(t,fluo3+4.5)
%     plot(t,fluo4+6)
%     plot(t,fluo5+7.5)
%     plot(t,fluo6+9)
%     plot(t,fluo6+10.5)
    
    plot(t,PeakArray_Speed-1.5,'k')
end


%%%%%%
filename_GCamp = ['dataMouseGCamp_',dataGCamp.Info.Name,'_',dataGCamp.Info.Date,'_Par'];
save(filename_GCamp,'dataGCamp')
%%%%%%




