%Script To compute Parameters from Force Peaks
%ComputeForceParScript

clear
close all
clc

%%%%  Choice of the animal and trial day
Animal_Name = 'GCaMP17_stroke_BoNT';
% Animal_Name = 'GCaMP18_stroke_BoNT';
%%%%

%%%%  Folder
UsbPort = 'I';
% AnimalMainDir = [UsbPort,':\LENS\Animals Data'];
% AnimalMainDir = ['C:\Users\Stefano\Google Drive\Piattaforma Stefano\LENS_SSSA\ELABORAZIONE DATA\_Per_Prove'];
AnimalMainDir = ['C:\Users\Stefano\Google Drive\Piattaforma Stefano\LENS_SSSA\ELABORAZIONE DATA\_data_MAT_GCamp_Store'];
AnimalCurrDir = [AnimalMainDir,'\',Animal_Name];
%%%%%%%%


%%%%%%%%
ListFolderAnimalDir = dir(AnimalCurrDir);


%for days
for lfcd=3:length(ListFolderAnimalDir)
    
    DayCurrDir           = ListFolderAnimalDir(lfcd,1).name;
    ListFolderDayCurrDir = dir([AnimalCurrDir,'\',DayCurrDir]);
    
%     filename = ['dataMouseGCamp_',Animal_Name,'_',DayCurrDir(1:2),'.mat'];
    filename = [DayCurrDir];
    
    if exist([AnimalCurrDir,'\',filename])
        
        %load gcamp
        load([AnimalCurrDir,'\',filename])
        
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
%                 for nFtr=1:length(i_F_Tr)
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
                        
                                          
%                 end
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
        
        %             ch = menu('Do you want to change threshold?','Yes','Not, it is ok');
        %             switch ch
        %                 case 1
        %                     whileSTOP = 0;
        %                     pause
        %
        %                 case 2
        %
        %                     %store in a new field of dataGCamp
        dataGCamp.PeaksForceNoCleanPeakPar = PeaksPar;
        dataGCamp.PeaksPar_Fx_Robot = [[1:NumTrials]', t_target, ForceMean, SubMov, Attempts, ForceMax];
        %                     whileSTOP = 1;
        %             end
        %
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% 3) Find start point of Force Peaks and Fluo Peaks from their first derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %num Fluorescence ROIs + Force Signal
    NumSignals = size(dataGCamp.fluoROI,2)+1;
    
    %cell to store info Peaks (first row -> force Peaks, second row -> corresponding Fluo Peaks in i-th ROI
    PeaksPar_Fx_Fluo = cell(2,NumSignals-1);
        
    %interval to consider
    IntFluoToPlot_sec   = [1 2]'; %[sec]
    IntFluoToPlot       = [IntFluoToPlot_sec(1) IntFluoToPlot_sec(2)]*Fs; % [points]
    
    %%% force  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %force peaks to take account of -> parto dai picchi di forza trovati
    ind_ForcePeaksToTake = sort(unique([find(dataGCamp.PeaksForceNoCleanPeakPar(:,Status_Curr_Par) == StatusOfInterest(1)); find(dataGCamp.PeaksForceNoCleanPeakPar(:,Status_Curr_Par) == StatusOfInterest(end))]));
    
    %force par (force is fixed) -> the peaks selected can change
    StartTime_ForcePar = dataGCamp.PeaksForceNoCleanPeakPar(ind_ForcePeaksToTake, StartTime_Curr_Par);
    %force duration
    Duration_ForcePar  = dataGCamp.PeaksForceNoCleanPeakPar(ind_ForcePeaksToTake, Duration_Curr_Par);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    id_old = 1;
    for id=2:NumSignals
        % for id=8:8
        
        %fluorescence sig
        FluoROI_Sig = dataGCamp.fluoROI(:,id-1);
        %filtering -> high Pass 0.5 Hz
        FluoROI_Sig = cheb2HPfilt(FluoROI_Sig,0.5,2,Fs);
        
        %for each peaks
        for iNPeaks=1:length(StartTime_ForcePar) %for each peaks
            
            %start Force
            stF = (StartTime_ForcePar(iNPeaks))*Fs;
            %end Force
            enF = (StartTime_ForcePar(iNPeaks)+Duration_ForcePar(iNPeaks))*Fs-1;
            
            %to align in agreement with 0 (t start of ForcePeak)
            stF_real = stF;
            
            
            %Force interval over threshold (no clean peak)
            IntervalForce_Sig = ForceSignal(stF_real:enF);
            %Force time (no clean peak)
            t_noCleanPeak     = [0:length(IntervalForce_Sig)-1]; % [points]
            t_noCleanPeak_sec = t_noCleanPeak/Fs; %[sec]
            
            
            %to check size of array (not to overcome border)
            dimOVS = (stF_real+ IntFluoToPlot(2)-length(ForceSignal)) * (stF_real+(IntFluoToPlot(2))-length(ForceSignal) > 0);
            
            %Force interval
            IntervalForce_Sig_Long   = ForceSignal(  stF_real-(IntFluoToPlot(1))    :   stF_real+(IntFluoToPlot(2))-dimOVS );
            
            %first derivative Force interval                      
            der_IntervalForce_Sig_Long = sgolayfilt(derivative(IntervalForce_Sig_Long,0.04),3,21);
            
            %Fluo interval
            IntervalFluoROI_Sig_Long = FluoROI_Sig(  stF_real-(IntFluoToPlot(1))    :   stF_real+(IntFluoToPlot(2))-dimOVS   );
            %first derivative Fluo interval
            %         der_IntervalFluoROI_Sig_Long = derivative(IntervalFluoROI_Sig_Long,0.04);
            der_IntervalFluoROI_Sig_Long = cheb2LPfilt(derivative(IntervalFluoROI_Sig_Long,0.04),5,3,Fs);
            
            %Fluo time
            t_long_sec = [t_noCleanPeak_sec(1) - IntFluoToPlot_sec(1)   : 1/Fs:   t_noCleanPeak_sec(1)+ IntFluoToPlot_sec(2)]'; %[sec]
            
            %interval pre-peak
            intPrePeak = [1:IntFluoToPlot(1)]';
            
            %threshold pre-onset Force (pre stF_real)
            ThPre_Force = median(der_IntervalForce_Sig_Long( intPrePeak))  + std(der_IntervalForce_Sig_Long( intPrePeak));
            %threshold pre-onset Fluo (pre stF_real)
            ThPre_Fluo  = median(der_IntervalFluoROI_Sig_Long( intPrePeak)) + std(der_IntervalFluoROI_Sig_Long( intPrePeak));
            
            
            %first derivative Force interval over ThPre_Force
            der_IntervalForce_Sig_Long_OT   = diff(der_IntervalForce_Sig_Long( intPrePeak)>ThPre_Force);
            StartForcePeak = find(der_IntervalForce_Sig_Long_OT==1)+1;
            
            if ~isempty(StartForcePeak)
                %real start of peak of force
                StartForcePeak                =  intPrePeak(StartForcePeak(end));   %[points]
                StartForcePeak_sec            = t_long_sec(StartForcePeak);         %[sec]
                BadForcePeak = 0;
            else
                display('bad Force peak (not starting from zero)')
                BadForcePeak = 1;
            end
            
            %first derivative Fluo interval over ThPre_Fluo
            der_IntervalFluoROI_Sig_Long_OT = diff(der_IntervalFluoROI_Sig_Long( intPrePeak)>ThPre_Fluo);
            StartFluoPeak = find(der_IntervalFluoROI_Sig_Long_OT==1)+1;
            
            if ~isempty(StartFluoPeak)
                %real start of peak of fluo
                StartFluoPeak                 =  intPrePeak(StartFluoPeak(end));   %[points]
                StartFluoPeak_sec             = t_long_sec(StartFluoPeak);         %[sec]
                BadFluoPeak = 0;
            else
                display('bad Force peak (not starting from zero)')
                BadFluoPeak = 1;
            end
            
            
            if ~(BadForcePeak) &&  ~(BadFluoPeak) %check if peaks are ok (it means that the peaks must have a positive first derivative - not flat (as example GCamp10_stroke_day5)
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% 2a) Find number of sub-peaks  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                          
                
                %%%% find positive peaks in Force derivative
                Int_Peak_Force    = [StartForcePeak:length(der_IntervalForce_Sig_Long)]';
                [pksValue pksLoc] = findpeaks(der_IntervalForce_Sig_Long(StartForcePeak:end));
                %
                %%%%%%%%%
                %%% Force peaks check
                %%%%%%%%
                prePoint = 5; preDiffValue = 0.1;
                i_delPk = [];
                for i=1:length(pksLoc)
                    pk    = (pksValue(i));
                    if (pksLoc(i)-prePoint)>0
                        prePoint_Ok = prePoint;
                    else
                        prePoint_Ok = 1;
                    end
                    
                    pkPre = (der_IntervalForce_Sig_Long(Int_Peak_Force(pksLoc(i)-prePoint_Ok)));
                    if (pk-pkPre)<=preDiffValue
                        i_delPk = [i_delPk; i];
                    end
                    
                end
                %%%%%%%%
                
                pksValue(i_delPk) = [];
                pksLoc(i_delPk) = [];
                %
                numPksDer_Force        = length(pksLoc);
                valuePksDer_Force      = pksValue;
                locPksDer_Force_sec    = t_long_sec(Int_Peak_Force(pksLoc)); %[sec]
                
                
                %%%% find positive peaks in Fluo derivative
                Int_Peak_Fluo    = [StartFluoPeak:length(der_IntervalFluoROI_Sig_Long)]';
                [pksValue pksLoc] = findpeaks(der_IntervalFluoROI_Sig_Long(StartFluoPeak:end));
                numPksDer_Fluo         = length(pksLoc);
                valuePksDer_Fluo       = pksValue;
                locPksDer_Fluo_sec     = t_long_sec(Int_Peak_Fluo(pksLoc)); %[sec]
                
                clear pksValue pksLoc
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %remove offset of FLUO
                IntervalFluoROI_Sig_Long = IntervalFluoROI_Sig_Long-min(IntervalFluoROI_Sig_Long);
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% 3) Compute Parameters         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %info peak
                trial_peak  = dataGCamp.PeaksForceNoCleanPeakPar(ind_ForcePeaksToTake(iNPeaks),Trial_Curr_Par);
                status_peak = dataGCamp.PeaksForceNoCleanPeakPar(ind_ForcePeaksToTake(iNPeaks),Status_Curr_Par);
                
                [ParForce ParFluo] = ParFromCorrespondingForceFluoPeaks(IntervalForce_Sig_Long, IntervalFluoROI_Sig_Long,...
                    StartForcePeak,         StartFluoPeak,...
                    numPksDer_Force,        numPksDer_Fluo,...
                    stF_real, IntFluoToPlot-dimOVS, trial_peak, status_peak, Fs, 1);
                
                if ParFluo(1,8)>35 %-> controlla la matrice MM (controllo sul valore dell'area sottesa al picco di Fluo (se troppo bassa probabilmenete non c'� picco di fluo corrispondente a picco force)
                    PeaksPar_Fx_Fluo{1,id-1} = [PeaksPar_Fx_Fluo{1,id-1}; ParForce];
                    PeaksPar_Fx_Fluo{2,id-1} = [PeaksPar_Fx_Fluo{2,id-1}; ParFluo];
                end
                
                
                
                %%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if 0
                    magfactor = 10;
                    figure
                    hold on
                    %Force long
                    plot(t_long_sec,IntervalForce_Sig_Long*magfactor,'b')
                    %Force no-clean peak
                    plot(t_noCleanPeak_sec,IntervalForce_Sig*magfactor,'g')
                    %derivative Force
                    plot(t_long_sec,der_IntervalForce_Sig_Long,'k')
                    %New Start Force Peak
                    scatter(StartForcePeak_sec,IntervalForce_Sig_Long(StartForcePeak)*magfactor,'k')
                    %Force Peaks
                    scatter(locPksDer_Force_sec,valuePksDer_Force,'k','Marker','s')
                    
                    %Fluo long
                    plot(t_long_sec,IntervalFluoROI_Sig_Long,'r')
                    %derivative Fluo
                    plot(t_long_sec,der_IntervalFluoROI_Sig_Long,'m')
                    %New Start Fluo Peak
                    scatter(StartFluoPeak_sec,IntervalFluoROI_Sig_Long(StartFluoPeak),'m')
                    %Fluo Peaks
                    scatter(locPksDer_Fluo_sec,valuePksDer_Fluo,'m','Marker','x')
                    
                    legend({'Force', 'ForceFirstPeak','deriv Force','New Start Force Peak', 'Force Peaks'...
                        'Fluo',                  'deriv Fluo', 'New Start Fluo Peak',  'Fluo Peaks'})
                    
                end
                
                
                if 0
                    %plot all the curves
                    
                    %first peak
                    if  iNPeaks==1 && id ~= id_old
                        figure('Name',['ROI n.',num2str(id-1)])
                        Med_to_Plot = [];
                    end
                    
                    %Fluo long
                    StartForcePeak_Whole_Point = ParForce(1,3)*Fs;
                    hold on
                    t_long_TO_PLOT_sec               = [- IntFluoToPlot_sec(1)   : 1/Fs:   IntFluoToPlot_sec(2)]'; %[sec]
                    IntervalFluoROI_Sig_Long_TO_PLOT = FluoROI_Sig(  (StartForcePeak_Whole_Point)-(IntFluoToPlot(1))    :   StartForcePeak_Whole_Point+(IntFluoToPlot(2))   );
                    plot(t_long_TO_PLOT_sec ,IntervalFluoROI_Sig_Long_TO_PLOT,'b')
                    id_old= id;
                    
                    Med_to_Plot = [Med_to_Plot IntervalFluoROI_Sig_Long_TO_PLOT];
                    
                    if iNPeaks == length(StartTime_ForcePar)
                        figure('Name',['Median in ROI n.',num2str(id-1)])
                        plot(t_long_TO_PLOT_sec ,median(Med_to_Plot,2),'b')
                        hold on
                        plot(t_long_TO_PLOT_sec ,median(Med_to_Plot,2)+std(Med_to_Plot,[],2),'r')
                        plot(t_long_TO_PLOT_sec ,median(Med_to_Plot,2)-std(Med_to_Plot,[],2),'r')
                        
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
            end% end if bad peak
            
            
        end %end for each peaks
        
        
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    dataGCamp.PeaksPar_Fx_Fluo  = PeaksPar_Fx_Fluo;
    display('Parameter computation end')
    
    
    %%% Save dataGCamp %%%
    %%%%%%
    dataGCamp.Info.Date = DayCurrDir(1:2);
    
    filename_GCamp = ['dataMouseGCamp_',dataGCamp.Info.Name,'_',dataGCamp.Info.Date,'_Par'];
%     filename_GCamp = [filename];
    save([AnimalCurrDir,'\',DayCurrDir,'\',filename_GCamp],'dataGCamp')
    
    
%     save(['C:\Users\CNR-SSSUP\Google Drive\Piattaforma Stefano\LENS_SSSA\ELABORAZIONE DATA\_data_MAT_GCamp_Store\',filename_GCamp],'dataGCamp')
    save(['C:\Users\Stefano\Google Drive\Piattaforma Stefano\LENS_SSSA\ELABORAZIONE DATA\_data_MAT_GCamp_Store\',filename_GCamp],'dataGCamp')
    
    %%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    else
        display([filename,' is missing'])
end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







