%Script To compute Parameters from Force Peaks
%ComputeForceParScript

% clear
% close all
clc

%%%%  Folder Info
UsbPort = 'I';
User = getenv('username');
%%%

%%%%  Choice of the animal and trial day
%%%%%%%%%%%
% ListAnimalTogether = {'GCaMPChR2_16_stroke_BoNT'};

RefOK = 1;
ListAnimalTogether = {'GCaMPChR2_7_control',  'GCaMPChR2_17_control', 'GCaMPChR2_18_control',...
                      'GCaMPChR2_20_control', 'GCaMPChR2_21_control', 'GCaMPChR2_23_control','GCaMPChR2_24_control',...
                      'GCaMPChR2_8_stroke',   'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke',...
                      'GCaMPChR2_3_stroke_BoNT',...
                      'GCaMP16_stroke_BoNT', 'GCaMP17_stroke_BoNT','GCaMP18_stroke_BoNT','GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT',...
                      'GCaMPChR2_13_stroke_BoNT',...
                      'GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT','GCaMPChR2_16_stroke_BoNT'};
                  
                  
                  
% RefOK = 0;                  
% ListAnimalTogether = {'GCaMP3_control', 'GCaMP4_control', 'GCaMPChR2_1_control',...
%                       'GCaMP9_stroke',  'GCaMP10_stroke', 'GCaMP11_stroke', 'GCaMP14_stroke','GCaMP15_stroke'};
%%%%%%%%%%%


%%%%
if RefOK == 1
    AnimalMainDir        = [UsbPort,':\LENS\Animals Data'];
    AnimalWhereSaveGCamp = ['C:\Users\',User,'\Google Drive\Piattaforma Stefano\LENS_SSSA\ELABORAZIONE DATA\_data_MAT_GCamp_Store'];
elseif RefOK == 0
    AnimalMainDir        = [UsbPort,':\LENS\Animals Data\NoBregmaREF'];
    AnimalWhereSaveGCamp = ['C:\Users\',User,'\Google Drive\Piattaforma Stefano\LENS_SSSA\ELABORAZIONE DATA\_data_MAT_GCamp_Store_NoBregmaREF'];
end
%%%%


for LATo = 1:length(ListAnimalTogether)
    
Animal_Name                    = [ListAnimalTogether{LATo}];

AnimalCurrDir                  = [AnimalMainDir,'\',Animal_Name];
AnimalWhereSaveGCampCurrDir    = [AnimalWhereSaveGCamp ,'\',Animal_Name];
sufx_name = '_Par_MIP';
%%%%%%%%


%%%%%%%%
ListFolderAnimalDir = dir(AnimalCurrDir);


%for days
for lfcd=3:length(ListFolderAnimalDir)
    
    DayCurrDir           = ListFolderAnimalDir(lfcd,1).name;
    display([Animal_Name,'_',DayCurrDir])
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
        %Speed signal
        SpeedSignal     =  dataGCamp.speed;
        
        
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
        IntFluoToPlot_sec   = [1 3]'; %[sec]
        IntFluoToPlot       = [IntFluoToPlot_sec(1) IntFluoToPlot_sec(2)]*Fs; % [points]
        
%         %%% force  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         %force peaks to take account of -> parto dai picchi di forza trovati
%         ind_ForcePeaksToTake = sort(unique([find(dataGCamp.PeaksForceNoCleanPeakPar(:,Status_Curr_Par) == StatusOfInterest(1)); find(dataGCamp.PeaksForceNoCleanPeakPar(:,Status_Curr_Par) == StatusOfInterest(end))]));
%         
%         %force par (force is fixed) -> the peaks selected can change
%         StartTime_ForcePar = dataGCamp.PeaksForceNoCleanPeakPar(ind_ForcePeaksToTake, StartTime_Curr_Par);
%         %force duration
%         Duration_ForcePar  = dataGCamp.PeaksForceNoCleanPeakPar(ind_ForcePeaksToTake, Duration_Curr_Par);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        id_old = 1;
        for id=1:NumSignals
            
            %fluorescence sig MIP
            FluoROI_Sig1 = dataGCamp.ROI_MIP.ROI_Signal(:,id);
            %filtering -> low pass Pass 9 Hz (to filter heart frequency)
            FluoROI_Sig = cheb2LPfilt(FluoROI_Sig1,9,2,Fs);
            
            
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
            
            
            %             FluoROI_Sig = cheb2HPfilt(FluoROI_Sig1,0.1,2,Fs);  %<--- NEW
            % %             FluoROI_Sig = cheb2LPfilt(FluoROI_Sig,10,2,Fs);   %<--- NEW
            
            %             FluoROI_Sig = cheb2HPfilt(FluoROI_Sig1,0.5,2,Fs);  %<--- NEW
            
            
            %             xFF = FluoROI_Sig;
            %             winDim = 4096;
            %             [ssAf f] = pwelch(xFF,winDim,winDim/2,[],Fs);
            %             figure
            %             plot(f,ssAf,'r');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            %for each peaks
            for iNPeaks=1:length(StartTime_ForcePar) %for each peaks
                
                iNPeaks = round(iNPeaks);
                
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
                der_IntervalForce_Sig_Long = derivative(IntervalForce_Sig_Long,0.04);
                %second derivative Force interval
                der2_IntervalForce_Sig_Long = derivative(der_IntervalForce_Sig_Long,0.04);
                
                %Fluo interval
                IntervalFluoROI_Sig_Long = FluoROI_Sig(  stF_real-(IntFluoToPlot(1))    :   stF_real+(IntFluoToPlot(2))-dimOVS   );
                %first derivative Fluo interval
                der_IntervalFluoROI_Sig_Long = derivative(IntervalFluoROI_Sig_Long,0.04);
                %second derivative Force interval
                der2_IntervalFluoROI_Sig_Long = derivative(der_IntervalFluoROI_Sig_Long,0.04);
                
                
                %Fluo time
                t_long     = [-(IntFluoToPlot(1))    :  (IntFluoToPlot(2))-dimOVS ]';
                t_long_sec = t_long/Fs;
%                 t_long_sec = [t_noCleanPeak_sec(1) - IntFluoToPlot_sec(1)   : 1/Fs:   t_noCleanPeak_sec(1)+ IntFluoToPlot_sec(2)]'; %[sec]
            

                %Speed interval
                IntervalSpeed_Sig_Long = SpeedSignal (  stF_real-(IntFluoToPlot(1))    :   stF_real+(IntFluoToPlot(2))-dimOVS   );
                
                
                %interval pre-peak
                intPrePeak = [1:IntFluoToPlot(1)]';
                
                
                %%%% ONSET OF THE PEAKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%% FORCE PEAK %%%%%%%%%
                
                % -> devo trovare dove la derivata prima è zero
                SeqBinaryDerForce       = der_IntervalForce_Sig_Long(intPrePeak)>0;
                SeqBinaryDerForce_Zero  = find(SeqBinaryDerForce==0);
                               
                %last value of zero in the array is the first zero crossing
                %of the first derivative = constant value
                StartForcePeak               = SeqBinaryDerForce_Zero(end);      %[points]
                StartForcePeak_sec           = t_long_sec(StartForcePeak);       %[sec]
                
                
                %%%% FLUO PEAK %%%%%%%%%
                
                SeqBinaryDerFLUO        = der_IntervalFluoROI_Sig_Long(intPrePeak)>0;
                SeqBinaryDerFLUO_Zero   = find(SeqBinaryDerFLUO==0);
                
                %nel caso il valore non sia proprio zero considero la derivata prima e vado a trovare il primo minimo dopo il valore massimo (== picco)) (che
                %corrisponde all'ultimo valore del vettore der_IntervalFluoROI_Sig_Long in cui ho cambiato segno e tolto l'offset del valore di picco (in modo tale da 
                %usare la fun findpeaks che cerca i massimi)
                if isempty(SeqBinaryDerFLUO_Zero)
                    [pksValue pksLoc] = findpeaks( -(der_IntervalFluoROI_Sig_Long(intPrePeak) - der_IntervalFluoROI_Sig_Long(intPrePeak(end))));
                    SeqBinaryDerFLUO_Zero = pksLoc;
                end
                
                
                StartFluoPeak                = SeqBinaryDerFLUO_Zero(end);       %[points]
                StartFluoPeak_sec            = t_long_sec(StartFluoPeak);        %[sec]
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %remove offset of FLUO
                IntervalFluoROI_Sig_Long = IntervalFluoROI_Sig_Long - IntervalFluoROI_Sig_Long(StartFluoPeak);
                %                     IntervalFluoROI_Sig_Long = IntervalFluoROI_Sig_Long-min(IntervalFluoROI_Sig_Long);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%% END OF THE PEAKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                
                %max of signal after onset but inside 2 sec
                InterPostStart          = StartForcePeak:StartForcePeak + 2*Fs-1;
                %to avoid overcoming time border
                InterPostStart(InterPostStart>length(IntervalFluoROI_Sig_Long)) = [];
                
                [maxFl imaxFl]          = max(IntervalFluoROI_Sig_Long(InterPostStart));
                %interval after max
                InterPostMax            = InterPostStart(imaxFl):length(IntervalFluoROI_Sig_Long);
                %binary signal after max (procedimento come sopra per lo start )
                SeqBinaryPostMax        = IntervalFluoROI_Sig_Long(InterPostMax)>0;
                SeqBinaryPostMax_Zero   = find(SeqBinaryPostMax==0);
                
                %nel caso il valore non sia proprio zero trovo il valore
                %minimo del segnale IntervalFluoROI_Sig_Long (non
                %posso usare la procedura di prima perchè posso avere più
                %picchi (e quindi minimi locali) nel profilo di della
                %derivata prima
                if length(InterPostMax)>3
                    
                    if isempty(SeqBinaryPostMax_Zero)
                        
                        [pksValue pksLoc] = findpeaks( -(IntervalFluoROI_Sig_Long(InterPostMax)- IntervalFluoROI_Sig_Long(InterPostMax(1))));
                        SeqBinaryPostMax_Zero = pksLoc;
                        
                        %se non trova nulla in 2 sec aumento la dimensione purché possa essere aumentata
                        if isempty(SeqBinaryPostMax_Zero)
                            
                            InterPostMax_long            = InterPostStart(imaxFl):IntFluoToPlot(2)-dimOVS;
                            if ~isempty(InterPostMax_long) && length(InterPostMax_long)>3
                                
                                [pksValue pksLoc]            = findpeaks( -(IntervalFluoROI_Sig_Long(InterPostMax_long)- IntervalFluoROI_Sig_Long(InterPostMax(1))));
                                
                                if isempty(SeqBinaryPostMax_Zero)
                                    SeqBinaryPostMax_Zero = 1;
                                    check_FluoA = 0;
                                else
                                    check_FluoA = 1;
                                end
                            else
                                SeqBinaryPostMax_Zero = 1;
                                check_FluoA = 0;
                            end
                            
                        else
                            check_FluoA = 1;
                        end
                        
                    else
                        check_FluoA = 1;
                    end
                    
                else
                    SeqBinaryPostMax_Zero = 1;
                    check_FluoA = 0;
                end
                
                EndFluoPeak             = (InterPostMax(SeqBinaryPostMax_Zero(1)));     %[points]
                EndFluoPeak_sec         = t_long_sec(EndFluoPeak);                      %[sec]
                
                
                if iNPeaks>1
                    %check%
                    %to avoid taking the same fluo peak more than one time
                    if (stF_PreviousPeak + Duration_FluoPeak_PreviousPeak < stF)
                        
                        check_FluoB = 1;
                    else
                        
                        check_FluoB = 0;
                        display('fluo peak already selected')
                    end
                else
                    
                    check_FluoB = 1;
                end
                
                
                Duration_FluoPeak_PreviousPeak = ( EndFluoPeak - StartFluoPeak + 1);  %[points]
                stF_PreviousPeak               = stF;
                
                
                %check sulla durata del picco di fluo
                if abs( EndFluoPeak -StartFluoPeak)<25 
                    
                    check_FluoC = 0;
                    display('fluo peak too short')
                else
                    check_FluoC = 1;
                end
                    
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%  check if fluo peaks have been already taken into account   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if check_FluoA && check_FluoB && check_FluoC 
                    
                    numPksDer_Force   = [];
                    numPksDer_Fluo    = [];
                         
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%    Compute Parameters         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %info peak
                    trial_peak  = dataGCamp.PeaksForceNoCleanPeakPar(ind_ForcePeaksToTake(iNPeaks),Trial_Curr_Par);
                    status_peak = dataGCamp.PeaksForceNoCleanPeakPar(ind_ForcePeaksToTake(iNPeaks),Status_Curr_Par);                  
                    
                    
                    [ParForce ParFluo] = ParFromCorrespondingForceFluoPeaks_v3(IntervalForce_Sig_Long, IntervalFluoROI_Sig_Long,...
                                                                               StartForcePeak,         StartFluoPeak,...
                                                                               [],                     EndFluoPeak,...
                                                                               numPksDer_Force,        numPksDer_Fluo,...
                                                                               stF_real, IntFluoToPlot,dimOVS, trial_peak, status_peak, Fs, 1);
                    
                    %store data
                    PeaksPar_Fx_Fluo{1,id} = [PeaksPar_Fx_Fluo{1,id}; ParForce];
                    PeaksPar_Fx_Fluo{2,id} = [PeaksPar_Fx_Fluo{2,id}; ParFluo];
                    
                   
                    
                else
                    display('fluo peak not taken')
                    
                end% end if bad peak
                
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
                    % %                         %second derivative Force
                    % %                         plot(t_long_sec,der2_IntervalForce_Sig_Long,'c')
                    %New Start Force Peak
                    scatter(StartForcePeak_sec,IntervalForce_Sig_Long(StartForcePeak)*magfactor,'k')
                    
                    %Fluo long
                    plot(t_long_sec,IntervalFluoROI_Sig_Long,'r')
                    %derivative Fluo
                    plot(t_long_sec,der_IntervalFluoROI_Sig_Long,'m')
                    % %                         %second derivative Force
                    % %                         plot(t_long_sec,der2_IntervalFluoROI_Sig_Long,'c')
                    %New Start Fluo Peak
                    scatter(StartFluoPeak_sec,IntervalFluoROI_Sig_Long(StartFluoPeak),'m')
                    %New End Fluo Peak
                    scatter(EndFluoPeak_sec,IntervalFluoROI_Sig_Long(EndFluoPeak),'r');
                    
                                        
                    legend({'Force', 'ForceFirstPeak','deriv Force','New Start Force Peak',...
                             'Fluo',                  'deriv Fluo', 'New Start Fluo Peak',})
                         
                    
                         
                    pos = dataGCamp.pos(stF_real-(IntFluoToPlot(1))    :   stF_real+(IntFluoToPlot(2))-dimOVS);
                    plot(t_long_sec,pos,'k')
                                        
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
                    t_long_TO_PLOT_sec               = [- IntFluoToPlot_sec(1)   : 1/Fs:   IntFluoToPlot_sec(2)-dimOVS/Fs]'; %[sec]
                    IntervalFluoROI_Sig_Long_TO_PLOT = FluoROI_Sig(  (StartForcePeak_Whole_Point)-(IntFluoToPlot(1))    :   StartForcePeak_Whole_Point+(IntFluoToPlot(2))-dimOVS   );
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

end %end LATo











