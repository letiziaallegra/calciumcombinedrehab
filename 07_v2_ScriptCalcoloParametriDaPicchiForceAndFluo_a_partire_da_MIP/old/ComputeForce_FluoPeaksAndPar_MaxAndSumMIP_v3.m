%Script To compute Parameters from Force Peaks
%ComputeForceParScript

clear
close all
clc

%%%%  Choice of the animal and trial day
%%%%%%%%%%%
Animal_Name = 'GCaMPChR2_7_control';
% Animal_Name = 'GCaMPChR2_17_control';
% Animal_Name = 'GCaMPChR2_18_control';
% Animal_Name = 'GCaMPChR2_20_control';
% Animal_Name = 'GCaMPChR2_21_control';

% Animal_Name = 'GCaMPChR2_8_stroke';
% Animal_Name = 'GCaMPChR2_9_stroke';
% Animal_Name = 'GCaMPChR2_19_stroke';
% Animal_Name = 'GCaMPChR2_22_stroke';

% Animal_Name = 'GCaMP16_stroke_BoNT';
% Animal_Name = 'GCaMP17_stroke_BoNT';
% Animal_Name = 'GCaMP18_stroke_BoNT';
% Animal_Name = 'GCaMPChR2_3_stroke_BoNT';
% Animal_Name = 'GCaMPChR2_11_stroke_BoNT';
% Animal_Name = 'GCaMPChR2_12_stroke_BoNT';
% Animal_Name = 'GCaMPChR2_13_stroke_BoNT';
% Animal_Name = 'GCaMPChR2_14_stroke_BoNT';
% Animal_Name = 'GCaMPChR2_15_stroke_BoNT';
% Animal_Name = 'GCaMPChR2_16_stroke_BoNT';
%%%%%%%%%%%
%%%%

%%%%  Folder
UsbPort = 'I';
AnimalMainDir        = [UsbPort,':\LENS\Animals Data'];
User = getenv('username');
AnimalWhereSaveGCamp = ['C:\Users\',User,'\Google Drive\Piattaforma Stefano\LENS_SSSA\ELABORAZIONE DATA\_data_MAT_GCamp_Store'];

AnimalCurrDir                  = [AnimalMainDir,'\',Animal_Name];
AnimalWhereSaveGCampCurrDir    = [AnimalWhereSaveGCamp ,'\',Animal_Name];
 sufx_name = '_Par_MIP';
%%%%%%%%


%%%%%%%%
ListFolderAnimalDir = dir(AnimalCurrDir);


%for days
for lfcd=3:length(ListFolderAnimalDir)
    
    DayCurrDir           = ListFolderAnimalDir(lfcd,1).name;
    ListFolderDayCurrDir = dir([AnimalCurrDir,'\',DayCurrDir]);
    
    DataFilename = ['dataMouseGCamp_',Animal_Name,'_',DayCurrDir(1:2),sufx_name,'.mat'];
    
    if exist([AnimalCurrDir,'\',DayCurrDir,'\',DataFilename])
        
        %load gcamp
        load([AnimalCurrDir,'\',DayCurrDir,'\',DataFilename])
        
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
        
%         %interval to consider
%         IntFluoToPlot_sec   = [1 2]'; %[sec]
%         IntFluoToPlot       = [IntFluoToPlot_sec(1) IntFluoToPlot_sec(2)]*Fs; % [points]
        
        
        %force signal
        ForceSignal     = -dataGCamp.fx;
        
        
        PeaksPar = dataGCamp.PeaksForceNoCleanPeakPar;
        
        %%% force  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %force peaks to take account of -> parto dai picchi di forza trovati
        ind_ForcePeaksToTake = sort(unique([find(dataGCamp.PeaksForceNoCleanPeakPar(:,Status_Curr_Par) == StatusOfInterest(1)); find(dataGCamp.PeaksForceNoCleanPeakPar(:,Status_Curr_Par) == StatusOfInterest(end))]));
        
        %force par (force is fixed) -> the peaks selected can change
        StartTime_ForcePar = dataGCamp.PeaksForceNoCleanPeakPar(ind_ForcePeaksToTake, StartTime_Curr_Par);
        %force duration
        Duration_ForcePar  = dataGCamp.PeaksForceNoCleanPeakPar(ind_ForcePeaksToTake, Duration_Curr_Par);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% 3) Find start point of Force Peaks and Fluo Peaks from their first derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %num Fluorescence ROIs
        NumSignals = size(dataGCamp.ROI_MIP.ROI_Signal,2);
        
        %cell to store info Peaks (first row -> force Peaks, second row -> corresponding Fluo Peaks in i-th ROI
        PeaksPar_Fx_Fluo = cell(2,NumSignals);
        
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
        for id=1:NumSignals
            
            %fluorescence sig MIP
            FluoROI_Sig = dataGCamp.ROI_MIP.ROI_Signal(:,id);
%             %filtering -> high Pass 0.5 Hz
%             FluoROI_Sig = cheb2HPfilt(FluoROI_Sig1,0.5,2,Fs);
 
            
%%%%%%%%%%%%%%%%%%%%%%%%%           
            %%%%notch
            %Design
% %             f_notch = 10.5;
% %             wo = f_notch/(Fs/2);
% %             bw = wo/35;
% %             [b,a] = iirnotch(wo,bw);
% %             %Application
% %             FluoROI_Sig2 = filtfilt(b, a, FluoROI_Sig);
%             xFF = FluoROI_Sig;
%             winDim = 4096;
%             [ssAf f] = pwelch(xFF,winDim,winDim/2,[],Fs);
%             figure
%             plot(f,ssAf);  


            FluoROI_Sig = cheb2HPfilt(FluoROI_Sig,0.1,2,Fs);  %<--- NEW
            FluoROI_Sig = cheb2LPfilt(FluoROI_Sig,15,2,Fs);   %<--- NEW
            
           
            
%             xFF = FluoROI_Sig;
%             winDim = 4096;
%             [ssAf f] = pwelch(xFF,winDim,winDim/2,[],Fs);
%             figure
%             plot(f,ssAf,'r');  

%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
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
                ThPre_Force_MedPre = median(der_IntervalForce_Sig_Long( intPrePeak));
                ThPre_Force        = ThPre_Force_MedPre  + std(der_IntervalForce_Sig_Long( intPrePeak));
                %threshold pre-onset Fluo (pre stF_real)
                ThPre_Fluo_MedPre = median(der_IntervalFluoROI_Sig_Long( intPrePeak));
                ThPre_Fluo  = ThPre_Fluo_MedPre + std(der_IntervalFluoROI_Sig_Long( intPrePeak));
                
                
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
                    display('bad Fluo peak (not starting from zero)')
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
                    IntervalFluoROI_Sig_Long = IntervalFluoROI_Sig_Long; % - ThPre_Fluo_MedPre; %min(IntervalFluoROI_Sig_Long);
                    
                    
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
                    
                    if ParFluo(1,8)>35 %-> controlla la matrice MM (controllo sul valore dell'area sottesa al picco di Fluo (se troppo bassa probabilmenete non c'è picco di fluo corrispondente a picco force)
                        PeaksPar_Fx_Fluo{1,id} = [PeaksPar_Fx_Fluo{1,id}; ParForce];
                        PeaksPar_Fx_Fluo{2,id} = [PeaksPar_Fx_Fluo{2,id}; ParFluo];
                    end
                    
                    
                    
                    %%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if 1
                        magfactor = 1;
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
                            figure('Name',['ROI n.',num2str(id)])
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
                            figure('Name',['Median in ROI n.',num2str(id)])
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
        
        
        dataGCamp.ROI_MIP_PeaksPar_Fx_Fluo  = PeaksPar_Fx_Fluo;
       
        
        %%% Save dataGCamp %%%
        %%%%%%
        DataFilename_Par = [DataFilename(1:end-4),'_Par'];
        save([AnimalCurrDir,'\',DayCurrDir,'\',DataFilename_Par],'dataGCamp')
        %%%%%%
        
        %%% Save in _dataGCamp -> backup per Analisi %%%
        save([AnimalWhereSaveGCampCurrDir,'\',DataFilename_Par],'dataGCamp')
        %%%
        
        display([DataFilename_Par,': Parameter computation end'])
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
        
    else
        display([DataFilename,' is missing'])
    end
    
end











