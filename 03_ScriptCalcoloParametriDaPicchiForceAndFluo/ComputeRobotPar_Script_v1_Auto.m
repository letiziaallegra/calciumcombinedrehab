%Script To compute Parameters from Force Peaks
%ComputeForceParScript

clear
close all
clc


%%%%  Folder
UsbPort = 'I';
% AnimalMainDir = [UsbPort,':\LENS\Animals Data'];
AnimalMainDir = [UsbPort,':\LENS\Animals Data\NoBregmaREF'];
% % AnimalMainDir = ['C:\Users\Stefano\Google Drive\Piattaforma Stefano\LENS_SSSA\ELABORAZIONE DATA\_Per_Prove'];
%%%%%%%%

%%%%%%%% Animal_Name %%%%%%%%%%%%%%%%%%%
% ListAnimalTogether = {  'GCaMPChR2_7_control', 'GCaMPChR2_17_control', 'GCaMPChR2_18_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
%                            'GCaMPChR2_20_control', 'GCaMPChR2_21_control',...
%                         'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke',...
%                         'GCaMP16_stroke_BoNT','GCaMP17_stroke_BoNT','GCaMP18_stroke_BoNT','GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT',...
%                            'GCaMPChR2_3_stroke_BoNT', 'GCaMPChR2_13_stroke_BoNT',... 
%                         'GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT','GCaMPChR2_16_stroke_BoNT'};
% 
ListAnimalTogether = {  'GCaMP3_control', 'GCaMP4_control',...
                        'GCaMP9_stroke','GCaMP10_stroke', 'GCaMP11_stroke', 'GCaMP14_stroke', 'GCaMP15_stroke',...
                        'GCaMPChR2_1_control'};

% %  
% LoadSubfix = '_Par';
% SaveSubfix = '_Par';
% % 
LoadSubfix = '_Par_MIP';
SaveSubfix = '_Par_MIP';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for LATo = 1:length(ListAnimalTogether)
    
    close all
    Animal_Name                    = [ListAnimalTogether{LATo}]
    
    %%%%%%%%
    AnimalCurrDir = [AnimalMainDir,'\',Animal_Name];
    ListFolderAnimalDir = dir(AnimalCurrDir);
    
    
    %for days
    for lfcd=3:length(ListFolderAnimalDir)
        
        DayCurrDir           = ListFolderAnimalDir(lfcd,1).name;
        ListFolderDayCurrDir = dir([AnimalCurrDir,'\',DayCurrDir]);
        
        filename = ['dataMouseGCamp_',Animal_Name,'_',DayCurrDir(1:2),LoadSubfix,'.mat'];
        %     filename = ['dataMouseGCamp_',Animal_Name,'_',DayCurrDir(1:2),'_Par_MIP_Par.mat'];
        
        if exist([AnimalCurrDir,'\',DayCurrDir,'\',filename])
            
            %load gcamp
            load([AnimalCurrDir,'\',DayCurrDir,'\',filename])
            
            %sampling frequency
            Fs = dataGCamp.Info.Fs;
            %task status
            Status = dataGCamp.status;
            %number of the trials (vector)
            TrialsNum = dataGCamp.InfoTrial.TrialsVector;
            %resting (no force/fluo peaks) interval -> used to quantify f0 for the normalization of fluo
            RestInterval = dataGCamp.Info.IntervalToFluoForceMean_PoinForceFreq;
            
            
            % status da considerare
            StatusOfInterest   = [3 4]';
            %columns
            Trial_Curr_Par     = 1;
            Status_Curr_Par    = 2;
            StartTime_Curr_Par = 3;
            Duration_Curr_Par  = 4;
            %
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%% 1) Find NoClean Force Peaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Force is now positive
            ForceSignal     = -dataGCamp.fx;
            Sig = ForceSignal;
            %Threshold 3std
            %%%% Threshold_Noise =  mean(RumoreCella(:,2))+ std(RumoreCella(:,2))*3
            Threshold_Noise = 0.0175;
            %%%% Threshold for peak force
            Threshold       = mean( abs(Sig(RestInterval(1):RestInterval(2)))) + 3*std(  abs(Sig(RestInterval(1):RestInterval(2)))  );
            
            %Speed
            Speed            = -dataGCamp.speed;
            Threshold_Speed  = mean(Speed)+std(Speed)*2;
            %
            
            %         %special case
            %         if strcmp(Animal_Name,'GCaMPChR2_8_stroke')
            %             Threshold_Noise = Threshold_Noise/10;
            %             Threshold       = Threshold/10;
            %             Threshold_Speed = 0.01; %GCaMPChR2_8_stroke
            %         end
            
            %find peaks over Threshold_Noise -> to compute Num of attempts
            [Inf_OverNoise] = PeaksFinder_v3(Sig ,Threshold_Noise, Fs,1);
            %find peaks over Threshold
            [Inf] = PeaksFinder_v3(Sig ,Threshold, Fs,1);
            
            %if peaks have been found over threshold -> to compute Num of attempts
            if ~isempty(Inf_OverNoise) %if peak over Noise (ON)
                PeakDurON = round(abs(Inf_OverNoise(:,2)*Fs));
                StartStopPeakON = [round(abs(Inf_OverNoise(:,1)*Fs))  round(abs(Inf_OverNoise(:,1)*Fs))+PeakDurON-1];
                % parameters from peaks
                [PeaksParON] = ParFromPeaks(ForceSignal, Fs,StartStopPeakON, TrialsNum, Status,1);
                
                %store in a new field of dataGCamp
                %             dataGCamp.PeakArray(:,iNS) = PeakArray;
                dataGCamp.PeaksForceNoCleanPeakPar_OverNoiseThreshold = PeaksParON;
                
            end
            
            %if peaks have been found
            if ~isempty(Inf) %if peak
                
                PeakDur = round(abs(Inf(:,2)*Fs));
                StartStopPeak = [round(abs(Inf(:,1)*Fs))  round(abs(Inf(:,1)*Fs))+PeakDur-1];
                
                PeakArray = zeros(length(Sig),1);
                for ai=1:size(StartStopPeak,1)
                    PeakArray(StartStopPeak(ai):StartStopPeak(ai)) = 1;
                end
                
                %
                %force rearrangement based on speed peaks
                %
                
                %speed
                [Inf_Speed] = PeaksFinder_v3(Speed,Threshold_Speed, Fs,1);
                %rearrange start and end peaks of speed
                [StartStopPeak_Speed PeakArray_Speed] = PeaksRearrange_v2(Speed,Inf_Speed,Fs,Speed,'off');
                
                if ~isempty(Inf_Speed) %if peak
                    [StartStopPeak_Force PeakArray_Force] = ForcePeaksCleaner_FromSpeed(StartStopPeak, StartStopPeak_Speed,Fs,PeakArray, Sig, Speed,'off');
                    StartStopPeak = StartStopPeak_Force;
                    PeakArray     = PeakArray_Force;
                else
                    error('no speed peaks')
                end
                
                
                %%%%% Find Info Peaks and Parameters %%%%%%%%%  %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%%
                
                % parameters from peaks
                [PeaksPar] = ParFromPeaks(ForceSignal, Fs,StartStopPeak, TrialsNum, Status,1);
                
                
                
            end %end if peak
            
            dataGCamp.PeaksForceNoCleanPeakPar = PeaksPar;
            dataGCamp.SpeedPeak                = Inf_Speed;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%% 2) Compute Parameters from the Robot   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            t_target  = [];
            ForceMean = [];
            SubMov    = [];
            Attempts  = [];
            ForceMax  = [];
            
            %%%
            t_target =  [dataGCamp.InfoTrial.Trials_Start_End_St3(:,2) - dataGCamp.InfoTrial.Trials_Start_End_St3(:,1)] /dataGCamp.Info.Fs;
            
            NumTrials = dataGCamp.InfoTrial.NumTrials;
            for ntf=1:NumTrials
                %find total num force peaks over threshold in that trial
                i_F_OT_Tr = find(dataGCamp.PeaksForceNoCleanPeakPar_OverNoiseThreshold(:,Trial_Curr_Par) == ntf);
                %find total num force peaks over threshold in that trial in status 3
                i_F_OT_Tr_St3  = sort(unique([find(dataGCamp.PeaksForceNoCleanPeakPar_OverNoiseThreshold(i_F_OT_Tr,Status_Curr_Par) == StatusOfInterest(1)); find(dataGCamp.PeaksForceNoCleanPeakPar_OverNoiseThreshold(i_F_OT_Tr,Status_Curr_Par) == StatusOfInterest(end))]));
                
                
                %total num force peaks
                NumTotForcePeaksOT = length(i_F_OT_Tr_St3);
                
                %find trial for peaks ok
                i_F_Tr = find(dataGCamp.PeaksForceNoCleanPeakPar(:,Trial_Curr_Par) == ntf);
                if ~isempty(i_F_Tr)
                    for nFtr=1:length(i_F_Tr)
                        %find status in trial for peaks ok
                        i_F_Tr_St3 = sort(unique([find(dataGCamp.PeaksForceNoCleanPeakPar(i_F_Tr,Status_Curr_Par) == StatusOfInterest(1)); find(dataGCamp.PeaksForceNoCleanPeakPar(i_F_Tr,Status_Curr_Par) == StatusOfInterest(end))]));
                        
                        %%%
                        ForceMean(ntf,1) =  mean(dataGCamp.PeaksForceNoCleanPeakPar(i_F_Tr(i_F_Tr_St3),5));
                        
                        %%%
                        SubMov(ntf,1)    =  length(dataGCamp.PeaksForceNoCleanPeakPar(i_F_Tr(i_F_Tr_St3),5));
                        
                        %%%
                        Attempts(ntf,1)  =  abs(NumTotForcePeaksOT-SubMov(ntf,1));
                        
                        %%%
                        FM  =  max(dataGCamp.PeaksForceNoCleanPeakPar(i_F_Tr(i_F_Tr_St3),5));
                        if ~isempty(FM)
                            ForceMax(ntf,1) = FM;
                        else
                            ForceMax(ntf,1) = NaN;
                        end
                        
                        
                    end
                else
                    %%%
                    ForceMean(ntf,1) =  NaN;
                    
                    %%%
                    SubMov(ntf,1)    =  NaN;
                    
                    %%%
                    Attempts(ntf,1)  =  NaN;
                    
                    %%%
                    ForceMax(ntf,1)  =  NaN;
                end
            end
            
            %figure Robot Parameters
            ParFig = figure;
            subplot(231)
            stem([1:NumTrials]',t_target);
            subplot(232)
            stem([1:NumTrials]',ForceMean);
            subplot(233)
            stem([1:NumTrials]',SubMov);
            subplot(234)
            stem([1:NumTrials]',Attempts);
            subplot(235)
            stem([1:NumTrials]',ForceMax);
            
            dataGCamp.PeaksForceNoCleanPeakPar = PeaksPar;
            dataGCamp.PeaksPar_Fx_Robot = [[1:NumTrials]', t_target, ForceMean, SubMov, Attempts, ForceMax];
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            display('Parameter computation end')
            
            
            %%% Save dataGCamp %%%
            %%%%%%
            
            filename_GCamp = ['dataMouseGCamp_',dataGCamp.Info.Name,'_',dataGCamp.Info.Date,SaveSubfix];
            %     filename_GCamp = [filename];
            save([AnimalCurrDir,'\',DayCurrDir,'\',filename_GCamp],'dataGCamp')
           
            %%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
        else
            display([filename,' is missing'])
        end
        
    end
    
end %end LATo


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







