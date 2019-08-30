%Script To compute Parameters from Force Peaks
%(using the already force peaks. These peaks correspond to visuo-selected Fluo sequences)
%ComputeForceParScript

%%%%%% last update %%%%%%
% 27/09/2016
%%%%%%%%%%%%%%%%%%%%%%%%%

%%% PREP CREATE DIRECTORY WITH ANIMAL NAME IN _data_MAT_GCamp_Store

%%% OUT A FILE PER DAY, MOVE THE RESULTS TO
%%% _data_MAT_GCamp_Store_Analysis_Par
clear
close all
clc

%%%%  Folder Info
UsbPort = 'M';
User = getenv('username');
%%%

%%%%  Choice of the animal and trial day
%%%%%%%%%%%
ListAnimalTogether = {'GCampChR2_TOX1','GCampChR2_TOX2','GCampChR2_TOX3','GCampChR2_TOX4','GCampChR2_TOX5'};
ListAnimalTogether = {  'GCamp_26_sani' };
ListAnimalTogether = {  'Tox1_toxin','Tox2_toxin','Tox3_toxin','Tox4_toxin','Tox5_toxin' }; 
ListAnimalTogether = {  'or19_robot','or20_robot','or21_robot','or23_optostim','or24_optostim',...
    'or25_optostim','or26_optostim+robot','or27_robot','or28_optostim+robot', 'or29_sham', 'or30_sham'};
ListAnimalTogether = {'GCaMP-ChR2-22_stroke', 'GCaMP20_Ctrl', 'GCaMP21_Ctrl' , 'GCaMP22_Ctrl', ...
    'GCaMP22_robot', 'GCaMP23_Ctrl', 'GCaMP23_robot', 'GCaMP24_robot', 'GCaMP25_robot', 'GCaMP26_robot', ...
    'GCaMP27_Ctrl', 'GCaMP28_Ctrl', 'GCaMP29_Ctrl'};
ListAnimalTogether = { 'GCaMP27_wk02_Ctrl', 'GCaMP27_wk03_Ctrl','GCaMP27_wk04_Ctrl','GCaMP28_wk02_Ctrl'...
,'GCaMP28_wk03_Ctrl','GCaMP28_wk04_Ctrl','GCaMP29_wk02_Ctrl','GCaMP29_wk03_Ctrl','GCaMP29_wk04_Ctrl'};
ListAnimalTogether = {'OR9_wk01_sham', 'OR10_wk01_sham', 'OR13_wk01_optostim','OR14_wk01_optostim'...
'OR15_wk01_optostim+robot','OR16_wk01_optostim+robot','OR17_wk01_optostim+robot'};
ListAnimalTogether = {'GCaMP24_wk01_ctrl','GCaMP25_wk01_ctrl','GCaMP26_wk01_ctrl'};
ListAnimalTogether = {'GCaMP-ChR2-1_wk01_ctrl'};

RefOK = 1;
% ListAnimalTogether = {'GCaMPChR2_7_control',  'GCaMPChR2_17_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
%                       'GCaMPChR2_8_stroke',   'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke', 'GCaMPChR2_25_stroke', 'GCaMPChR2_26_stroke',...
%                       'GCaMP16_stroke_BoNT','GCaMP18_stroke_BoNT','GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT',...
%                       'GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT','GCaMPChR2_16_stroke_BoNT'};

% ListAnimalTogether = {'GCaMPChR2_32_stroke_Rehab'};

% %
% RefOK = 0;
% ListAnimalTogether = {'GCaMP3_control', 'GCaMP4_control', 'GCaMPChR2_1_control',...
%                       'GCaMP9_stroke',  'GCaMP10_stroke', 'GCaMP11_stroke', 'GCaMP14_stroke','GCaMP15_stroke'};
%%%%%%%%%%%


%%%%
if RefOK == 1
    %AnimalMainDir       = [UsbPort,':\LENS\Animals Data STIM'];
    AnimalMainDir        = ['C:\LENS\Data Leti'];
    AnimalMainDir        = ['/Volumes/ALE6TB_DESK/DATA/Imaging/Emilia/Rehab/MATLAB_DATA_FOLDERS/'];
    AnimalMainDir        = '/Users/alessandro/Desktop/toxin/MATLAB/';
    AnimalMainDir        = '/Users/alessandro/Desktop/180424_RehabOptogen/MATLAB/';
    AnimalMainDir        = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB/';
    AnimalMainDir        = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_GCaMP27-29/';
    AnimalMainDir        = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_OR';
    AnimalMainDir        = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_CTRL_190122';
    AnimalMainDir        = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_CTRL_190117';
    
    
    AnimalWhereSaveGCamp = ['C:\Users\asus\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\_data_MAT_GCamp_Store'];
    AnimalWhereSaveGCamp = ['/Users/alessandro/Desktop/ELABORAZIONE DATA/_data_MAT_GCamp_Store'];
elseif RefOK == 0
    AnimalMainDir        = [UsbPort,':\LENS\Animals Data\NoBregmaREF'];
    AnimalWhereSaveGCamp = ['C:\Users\',User,'\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\_data_MAT_GCamp_Store'];
end
%%%%


for LATo = 1:length(ListAnimalTogether)
    
    Animal_Name                    = [ListAnimalTogether{LATo}];
    
    AnimalCurrDir                  = [AnimalMainDir,filesep,Animal_Name];
    AnimalWhereSaveGCampCurrDir    = [AnimalWhereSaveGCamp ,filesep,Animal_Name];
    sufx_name = '_Par_MIPSIP';
    %%%%%%%%
    
    
    %%%%%%%%
    ListFolderAnimalDir = dir(AnimalCurrDir);
    
    
    %for days
    for lfcd=3:length(ListFolderAnimalDir)
        
        DayCurrDir           = ListFolderAnimalDir(lfcd,1).name;
        display([Animal_Name,'_',DayCurrDir])
        ListFolderDayCurrDir = dir([AnimalCurrDir,filesep,DayCurrDir]);
        
        DataFilename = ['dataMouseGCamp_',Animal_Name,'_',DayCurrDir(1:2),sufx_name,'.mat'];
        
        if exist([AnimalCurrDir,filesep,DayCurrDir,filesep,DataFilename])
            
            %load GCamp file
            load([AnimalCurrDir,filesep,DayCurrDir,filesep,DataFilename])
            %load file info Seq -> DataInfoSequence
            load([AnimalCurrDir,filesep,DayCurrDir,filesep,'SequenceLong',filesep,'DataInfoSequence_',Animal_Name,'_',DayCurrDir(1:2)])
            SequeOK = find(DataInfoSequence(:,3)==1);
            
            
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
            %Pos signal
            PosSignal       =  dataGCamp.pos;
            %Speed signal
            SpeedSignal     =  dataGCamp.speed;
            %Acceleration signal
            AccSignal       =  dataGCamp.acc;
            
            
            
            %%% force  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%% N.B. the "good" force peaks have been already selected  %%%
            SelectedForcePeaks = dataGCamp.PeaksPar_Fx_Fluo{1,1}; % -> force peaks found in the previous analysis (whole image)
            GoodForcePeaks     = SelectedForcePeaks(SequeOK,:);   % -> among these force peaks found in the previous analysis (whole image), we visually selected a fraction of them
            
            
            %         %force peaks to take account of -> parto dai picchi di forza trovati
            %         ind_ForcePeaksToTake = sort(unique([find(dataGCamp.PeaksForceNoCleanPeakPar(:,Status_Curr_Par) == StatusOfInterest(1)); find(dataGCamp.PeaksForceNoCleanPeakPar(:,Status_Curr_Par) == StatusOfInterest(end))]));
            %
            %onset "good" force peaks
            StartTime_ForcePar = GoodForcePeaks(:, StartTime_Curr_Par);
            %duration "good" force peaks
            Duration_ForcePar  = GoodForcePeaks(:, Duration_Curr_Par);
            %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%% 3) Find start point of   ONLY  Fluo Peaks         from their first derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %num Fluorescence ROIs
            NumSignals = size(dataGCamp.ROI_MIP_SIP.ROI_Signal,2);
            
            %cell to store info Peaks (first row -> force Peaks, second row -> corresponding Fluo Peaks in i-th ROI
            PeaksPar_Fx_Fluo = cell(2,NumSignals);
            %cell to store fluo curve aligned to start movement
            FluoPeaksAligned_to_FirstMov = cell(1,NumSignals);
            
            %interval to consider
            IntFluoToPlot_sec   = [1 3]'; %[sec]
            IntFluoToPlot       = [IntFluoToPlot_sec(1) IntFluoToPlot_sec(2)]*Fs; % [points]
            
            id_old        = 1;
            t_to_plot_all = [-IntFluoToPlot(1) : IntFluoToPlot(2)]';
            
            
            %%%%% Scroll ROI-FLUO Signals %%%%%%%%
            for id=1:NumSignals
                
                %fluorescence sig MIP or SIP
                FluoROI_Sig1 = dataGCamp.ROI_MIP_SIP.ROI_Signal(:,id);
                
                %filtering -> low pass Pass 9 Hz (to filter heart frequency)
                FluoROI_Sig  = cheb2LPfilt(FluoROI_Sig1,9,2,Fs);
                
                
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
                iNPeaks_all = 0;
                for iNPeaks=1:size(GoodForcePeaks,1) %for each force peaks
                    
                    iNPeaks = round(iNPeaks);
                    
                    
                    %%%% DATA  PEAKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    %start Force
                    stF = round((StartTime_ForcePar(iNPeaks))*Fs);
                    %end Force
                    enF = (StartTime_ForcePar(iNPeaks)+Duration_ForcePar(iNPeaks))*Fs-1;
                    
                    %to align in agreement with 0 (t start of ForcePeak)
                    stF_real = stF;
                    
                    %to align in agreement with the interval (it is the point in the middle)
                    StartForcePeak_REL_INT = IntFluoToPlot_sec(1)+1 *Fs; %[points]  -> related to the interval (it is central point of the interval (the IntFluoToPlot_sec+1 point))
                    enF = round(enF);
                    %Force interval over threshold (no clean peak)
                    IntervalForce_Sig = ForceSignal(stF_real:enF);
                    %Force time (no clean peak)
                    t_noCleanPeak     = [0:length(IntervalForce_Sig)-1]; % [points]
                    t_noCleanPeak_sec = t_noCleanPeak/Fs; %[sec]
                    
                    
                    %to check size of array (not to overcome border)
                    dimOVS = (stF_real+ IntFluoToPlot(2)-length(ForceSignal)) * (stF_real+(IntFluoToPlot(2))-length(ForceSignal) > 0);
                    
                    
                    %Force interval
                    IntervalForce_Sig_Long      = ForceSignal(  stF_real-(IntFluoToPlot(1))    :   stF_real+(IntFluoToPlot(2))-dimOVS );
                    %first derivative Force interval
                    der_IntervalForce_Sig_Long  = derivative(IntervalForce_Sig_Long,0.04);
                    %second derivative Force interval
                    der2_IntervalForce_Sig_Long = derivative(der_IntervalForce_Sig_Long,0.04);
                    
                    
                    %Fluo interval
                    IntervalFluoROI_Sig_Long = FluoROI_Sig(  stF_real-(IntFluoToPlot(1))    :   stF_real+(IntFluoToPlot(2))-dimOVS   );
                    %first derivative Fluo interval
                    der_IntervalFluoROI_Sig_Long = derivative(IntervalFluoROI_Sig_Long,0.04);
                    %second derivative Force interval
                    der2_IntervalFluoROI_Sig_Long = derivative(der_IntervalFluoROI_Sig_Long,0.04);
                    
                    
                    %Pos interval
                    IntervalPos_Sig_Long   = PosSignal   (  stF_real-(IntFluoToPlot(1))    :   stF_real+(IntFluoToPlot(2))-dimOVS   );
                    %Speed interval
                    IntervalSpeed_Sig_Long = SpeedSignal (  stF_real-(IntFluoToPlot(1))    :   stF_real+(IntFluoToPlot(2))-dimOVS   );
                    %Speed interval
                    IntervalAcc_Sig_Long   = AccSignal   (  stF_real-(IntFluoToPlot(1))    :   stF_real+(IntFluoToPlot(2))-dimOVS   );
                    
                    
                    %Fluo time
                    t_long     = [-(IntFluoToPlot(1))    :  (IntFluoToPlot(2))-dimOVS ]';
                    t_long_sec = t_long/Fs;
                    %                 t_long_sec = [t_noCleanPeak_sec(1) - IntFluoToPlot_sec(1)   : 1/Fs:   t_noCleanPeak_sec(1)+ IntFluoToPlot_sec(2)]'; %[sec]
                    
                    %interval pre-peak
                    intPrePeak = [1:IntFluoToPlot(1)]';
                    
                    
                    %%%% ONSET OF THE PEAKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    %%%% FORCE PEAK %%%%%%%%%
                    
                    %already found
                    StartForcePeak               = StartForcePeak_REL_INT;      %[points] -> related to the whole signal
                    StartForcePeak_sec           = (StartForcePeak-IntFluoToPlot(1))/Fs;           %[sec]    -> related to the whole signal
                    
                    
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
                    
                    StartFluoPeak                = SeqBinaryDerFLUO_Zero(end);       %[points]  -> related to the interval
                    StartFluoPeak_sec            = t_long_sec(StartFluoPeak);        %[sec]     -> related to the interval
                    
                    
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
                    
                    check_FluoC = 1;
                    
                    
                    
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
                        trial_peak  = GoodForcePeaks(iNPeaks,Trial_Curr_Par);
                        status_peak = GoodForcePeaks(iNPeaks,Status_Curr_Par);
                        
                        if (Duration_FluoPeak_PreviousPeak)>3
                            [ParForce ParFluo] = ParFromCorrespondingForceFluoPeaks_v4(IntervalForce_Sig_Long, IntervalFluoROI_Sig_Long,  IntervalAcc_Sig_Long,...
                                StartForcePeak,         StartFluoPeak,...
                                [],                     EndFluoPeak,...
                                numPksDer_Force,        numPksDer_Fluo,...
                                stF_real, IntFluoToPlot,dimOVS, trial_peak, status_peak, Fs, 1);
                            
                            %store data
                            PeaksPar_Fx_Fluo{1,id} = [PeaksPar_Fx_Fluo{1,id}; ParForce];
                            PeaksPar_Fx_Fluo{2,id} = [PeaksPar_Fx_Fluo{2,id}; ParFluo];
                        else
                            display('interval too short')
                        end
                        
                        
                    else
                        display('fluo peak not taken')
                        
                    end% end if bad peak
                    
                    %%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if check_FluoA && check_FluoB && check_FluoC
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
                            
                            %Fluo long
                            plot(t_long_sec,IntervalFluoROI_Sig_Long,'r')
                            %derivative Fluo
                            plot(t_long_sec,der_IntervalFluoROI_Sig_Long,'m')
                            %New Start Fluo Peak
                            scatter(StartFluoPeak_sec,IntervalFluoROI_Sig_Long(StartFluoPeak),'m')
                            %New End Fluo Peak
                            scatter(EndFluoPeak_sec,IntervalFluoROI_Sig_Long(EndFluoPeak),'r');
                            
                            %Pos
                            plot(t_long_sec,IntervalPos_Sig_Long/10,'b')
                            %Acc
                            plot(t_long_sec,IntervalAcc_Sig_Long/50,'c')
                            %First Mov
                            istFirstMov = round( (ParForce(1,13)+ ParForce(1,3))*Fs - stF_real + IntFluoToPlot(1) );
                            scatter(t_long_sec(istFirstMov),IntervalPos_Sig_Long(istFirstMov)/10,'r')
                            
                            legend({'Force', 'ForceFirstPeak','deriv Force','New Start Force Peak',...
                                'Fluo',                  'deriv Fluo' , 'New Start Fluo Peak',  'End Fluo Peak', ...
                                'Pos/10','Acceleration/50'})
                            
                        end
                        
                        
                        if 1
                            %Med_to_Plot
                            iNPeaks_all = iNPeaks_all+1;
                            
                            %first peak
                            if  iNPeaks_all ==1
                                Med_to_Plot = [];
                            end
                            
                            %start position
                            Start_Whole_Point = round((ParForce(1,13)+ ParForce(1,3))*Fs - stF_real + IntFluoToPlot(1));
                            
                            %Fluo sig
                            Fluo_to_plot_all      = IntervalFluoROI_Sig_Long(StartFluoPeak:EndFluoPeak);
                            %time Fluo Sig
                            t_Fluo_to_plot_all    = [StartFluoPeak-Start_Whole_Point:EndFluoPeak-Start_Whole_Point]';
                            
                            
                            %find time to align
                            ind_T = find(t_to_plot_all == t_Fluo_to_plot_all(1));
                            
                            if ~isempty(ind_T)
                                
                                %matrix the current peak
                                Med_to_Plot_buf  = ones(length(t_to_plot_all),1) *NaN;
                                
                                ind_T_end = ind_T+length(t_Fluo_to_plot_all)-1;
                                max_T_end = IntFluoToPlot_sec(1)*Fs + 1 + IntFluoToPlot_sec(2)*Fs -dimOVS;
                                if ind_T_end  > max_T_end
                                    dif_len = (ind_T_end-max_T_end);
                                    ind_T_end = ind_T_end - dif_len;
                                else
                                    dif_len = 0;
                                end
                                
                                
                                
                                Med_to_Plot_buf(ind_T:ind_T_end,1) = Fluo_to_plot_all(1:end-dif_len);
                                
                                %matrix with all the peaks
                                Med_to_Plot = [Med_to_Plot, Med_to_Plot_buf];
                                
                            end
                        end
                        
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                end %end for each peaks
                
                %%%%%%%%%%%%%%%%Save Med_to_Plot
                FluoPeaksAligned_to_FirstMov{1,id} = Med_to_Plot;
                if 0
                    figure('Name',['Median in ROI n.',num2str(id)])
                    plot(t_to_plot_all/Fs ,nanmedian(Med_to_Plot,2),'b')
                    hold on
                    plot(t_to_plot_all/Fs ,nanmedian(Med_to_Plot,2)+nanstd(Med_to_Plot,[],2),'r')
                    plot(t_to_plot_all/Fs ,nanmedian(Med_to_Plot,2)-nanstd(Med_to_Plot,[],2),'r')
                end
                %%%%%%%%%%%%%%%%
                
                
            end %end ROI
            close all
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            dataGCamp.ROI_MIP_SIP_PeaksPar_Fx_Fluo              = PeaksPar_Fx_Fluo;
            dataGCamp.ROI_MIP_SIP_FluoPeaksAligned_to_FirstMov  = FluoPeaksAligned_to_FirstMov;
            
            
            %%% Save dataGCamp %%%
            if 1
                
                DataFilename_Par = [DataFilename(1:end-4),'_Par'];
                save([AnimalCurrDir,filesep,DayCurrDir,filesep,DataFilename_Par],'dataGCamp')
                %%%%%%
                
                %%% Save in _dataGCamp -> backup per Analisi %%%
                save([AnimalWhereSaveGCampCurrDir,filesep,DataFilename_Par],'dataGCamp')
                %%%
                
                display([DataFilename_Par,': Parameter computation end'])
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
        else
            display([DataFilename,' is missing'])
        end
        
    end
    
end %end LATo

display('Process completed')









