%
% script to syncrhonize fluo and force signals and
% define the Fluo ROIs signals
%
clear
close all
clc

%% Choice of the animal and trial day
%UsbPort = 'M';
%AnimalDir = [UsbPort,':\LENS\Animals Data STIM'];
AnimalDir = ['C:\LENS\Data Leti'];
% AnimalDir = ['C:\Users\CNR-SSSUP\Desktop\DATI LENS\Animals_Data_Make_Mat_Seq_to_do'];

%%%%%%%% Animal_Name
% Animal_Name_choice = 'GCaMPChR2_7_control';
% Animal_Name_choice = 'GCaMPChR2_8_stroke';
% Animal_Name_choice = 'GCaMPChR2_3_stroke_BoNT';
% Animal_Name_choice = 'GCaMP16_stroke_BoNT';
% Animal_Name_choice = 'GCaMP17_stroke_BoNT';
% Animal_Name_choice = 'GCaMP18_stroke_BoNT';
% Animal_Name_choice = 'GCaMPChR2_14_stroke_BoNT';
% Animal_Name_choice = 'GCaMPChR2_24_control';
% Animal_Name_choice = 'GCaMPChR2_26_stroke_LateRehab'; 

Animal_Name_choice = 'GCampChR2_TOX1';

%%%%%%%% TrialDay
% TrialDay_choice = '01';
TrialDay_choice = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%% choice of the animals to elaborate     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ListFolderAnimalDir = dir(AnimalDir);
indexAn = 1;
for lfcd=3:length(ListFolderAnimalDir)
    AnimalName_buf = ListFolderAnimalDir(lfcd,1).name;
    
    if length(AnimalName_buf)>5
        
        if strcmp(AnimalName_buf(1:5),'GCaMP') & isdir([AnimalDir,'\',AnimalName_buf])
            %all of the animals in folder
            AnimalName_IndexList(indexAn,1) = lfcd;
            indexAn = indexAn+1;
            
        end
        
        if strcmp(AnimalName_buf,Animal_Name_choice) & isdir([AnimalDir,'\',AnimalName_buf])
            %one animal
            AnimalName_IndexList_One = lfcd;
            break
        end
    end
end
if ~isempty(Animal_Name_choice)
    %one animal
    AnimalName_IndexList = AnimalName_IndexList_One;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for anAn_i=1:length(AnimalName_IndexList)  %for Animale -> end
        
    Animal_Index = AnimalName_IndexList(anAn_i);
    Animal_Name  = ListFolderAnimalDir(Animal_Index,1).name;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% choice of the days to elaborate        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AnimalCurrDir = dir([AnimalDir,'\',Animal_Name]);
    indexDay = 1;
    
    for acd=3:length(AnimalCurrDir)
        AnimalDay_buf = AnimalCurrDir(acd,1).name;
        
        if isdir([AnimalDir,'\',Animal_Name,'\',AnimalDay_buf])
            %all of the days of the animals in folder
            Days_IndexList(indexDay,1) = acd;
            indexDay = indexDay+1;
            
            if strcmp(AnimalDay_buf(1:2),TrialDay_choice)
                %one animal
                Days_IndexList_One = acd;
                break
            end
        end
    end
    if ~isempty(TrialDay_choice)
        %one animal
        Days_IndexList = Days_IndexList_One;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %scroll days
    for anD_i=1:length(Days_IndexList)
        
        TrialDay_Index = Days_IndexList(anD_i);
        TrialDay       = AnimalCurrDir(TrialDay_Index,1).name;
        
        CurrAnDayFolder = [AnimalDir,'\',Animal_Name,'\',TrialDay];
        CurrAnDayFolder_List = dir(CurrAnDayFolder);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% find image folder and force file       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        folderTASK_FLUO          = [];
        folderTASK_ForceFileName = [];
        for cadf=3:length(CurrAnDayFolder_List)
            
            if length(CurrAnDayFolder_List(cadf,1).name)>2
                
                if strcmp(CurrAnDayFolder_List(cadf,1).name(1:3),'MAT')
                    %data images folder
                    folderTASK_FLUO = [AnimalDir,'\',Animal_Name,'\',TrialDay,'\',CurrAnDayFolder_List(cadf,1).name];
                elseif length(CurrAnDayFolder_List(cadf,1).name)>7
                    if strcmp(CurrAnDayFolder_List(cadf,1).name(end-7:end-4),'sync')
                        %force file name
                        folderTASK_ForceFileName = [AnimalDir,'\',Animal_Name,'\',TrialDay,'\',CurrAnDayFolder_List(cadf,1).name];
                    end
                end
                
            end
            
        end
        if isempty(folderTASK_FLUO) && isempty(folderTASK_ForceFileName)
            error([folderTASK_FLUO,' or ', folderTASK_ForceFileName, 'missing']);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %size ROIs
        ROI_name = {'Whole'};
        %region of interest to take
        ROI_num  = length(ROI_name);
        %%
        
        %load force file
        dataRobot = importdata(folderTASK_ForceFileName);
        
        
        %%%%%% data Robot %%%%%%%%%%%%
        
        synchroTask = find(dataRobot(:,7)==0);
        %To synchronize robot file (the previous time point)
        synchroTask = synchroTask(1)-1; %-215;
        
        %starting point
        RealStart = 1;
        RealStart_Robot = synchroTask+RealStart;
        %duration point
        RealDurTask_Robot = size(dataRobot,1)-RealStart_Robot+1;
        % RealDurTask_Robot = 15000; %-> in 150729_GCaMP9_trailx15_8bits sono circa 4 trials
        fs_robot = 100; %Hz
        
        t_robot = dataRobot(RealStart_Robot:RealStart_Robot+RealDurTask_Robot-1,1);
        t_robot = t_robot - t_robot(1);
        fx_robot = dataRobot(RealStart_Robot:RealStart_Robot+RealDurTask_Robot-1,2);
        %filtering
        fx_robot = sgolayfilt(fx_robot-median(fx_robot),3,21);
        %status
        status_robot = dataRobot(RealStart_Robot:RealStart_Robot+RealDurTask_Robot-1,7);
        %
        pos_robot = [];
        if size(dataRobot,2)>8
            %position
            pos_robot = dataRobot(RealStart_Robot:RealStart_Robot+RealDurTask_Robot-1,9);
            %speed
            speed_robot = dataRobot(RealStart_Robot:RealStart_Robot+RealDurTask_Robot-1,10);
            %acceleration
            acc_robot = dataRobot(RealStart_Robot:RealStart_Robot+RealDurTask_Robot-1,11);
        end
        %%%%%%
        
        %%%%%% Take task images  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        d = dir(folderTASK_FLUO);
        durationinframes=10;
        len = length(d);
        NumImages = len-2;
        fs_fluo = 25; %Hz
        
        %resampling fsample Robot Fluo
        ResampleParRobotFluo = round(fs_robot/fs_fluo);
        %starting point
        RealStart_Fluo=RealStart;
        %duration
        RealDurTask_Fluo = floor(RealDurTask_Robot/ResampleParRobotFluo);
        
        
        %%%%%  check length of the files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if NumImages<RealDurTask_Fluo
            %difference
            diflen = RealDurTask_Fluo - NumImages;
            %differenza espressa in numero di punti del file del robot
            diflenR = diflen*ResampleParRobotFluo;
            RealDurTask_Fluo = NumImages;
        else
            diflenR = 0;
        end
        
        %risistemo files
        modToremov   = mod(RealDurTask_Robot,ResampleParRobotFluo);
        
        t_robot      = t_robot(1:end-diflenR-modToremov);
        fx_robot     = fx_robot(1:end-diflenR-modToremov);
        status_robot = status_robot(1:end-diflenR-modToremov);
        
        if ~isempty(pos_robot)
            pos_robot   = pos_robot(1:end-diflenR-modToremov);
            speed_robot = speed_robot(1:end-diflenR-modToremov);
            acc_robot   = acc_robot(1:end-diflenR-modToremov);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        downsamplingfactor=1;
        
        index=0;
        wb = waitbar(0,'Images Loading, Please wait...');
        
        % % % % %load image stack
        % % % % % N.B. size of StoreImage depends on num of ROIs
        % % % % StoreImageAllROIs = cell(ROI_num,1);
        % % % % ROI_PosAndSize    = zeros(ROI_num,4);
        
        
        
        LenFLUO_Im = RealStart_Fluo+RealDurTask_Fluo-1;
        for i=RealStart_Fluo :downsamplingfactor:   LenFLUO_Im
            
            index = index+1;
            indexImage = i+2;
%             waitbar(indexImage/(LenFLUO_Im),wb)
            display( num2str(indexImage/(LenFLUO_Im)))
            
            nameImage = d(indexImage,1).name;

            load([folderTASK_FLUO,'\',nameImage]);
            Im = Im8_fv;
            Im_Original = Im;          
            
            %filtering
%             Im = medfilt2(Im);
            StoreImageAllROIs(:,:,index) = Im;
            
        end
        
        %save size of the images
        ROI_PosAndSize = [1 1 size(Im,1)-1 size(Im,2)-1];
        
        clear d dataRobot Im Im_Original
        close(wb)
        
        %%%% choose the baseline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%
        
        
        %%%% plot time course  %%%
        fbas = figure('Name','Select where animal is not moving');
        plot(t_robot,status_robot/10,'r')
        hold on
        plot(t_robot,fx_robot,'b')
        
        StoreMean    = zeros(LenFLUO_Im);
        
        %compute average on the areas
        StoreMean    = squeeze(mean(mean(StoreImageAllROIs(:,:,:),2),1));
        ImFLUO_buf   = resample(StoreMean,ResampleParRobotFluo,1);
        
        %rimuove drift        
        plot(t_robot,(ImFLUO_buf - median(ImFLUO_buf))/10,'g');
        legend([{'Status','Force'},ROI_name])
        ylim([-3 3])
        
        [x_1 x_2] = ginput;
        %estremi intervallo di interesse
        IntervalToFluoMean           = round([x_1(1) x_1(2)]*fs_fluo);
        IntervalToFluoMean_ForceFreq = round([x_1(1) x_1(2)]);
        close(fbas)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        StoreImageMean = round(mean(double(StoreImageAllROIs(:,:,IntervalToFluoMean(1): IntervalToFluoMean(2))),3));
        StoreImageMean(StoreImageMean==0) = 1;
        
        %normalization
        for j=1:size(StoreImageAllROIs,3)
            %(f-f0)/f0
            StoreImageAllROIs(:,:,j) = ((double(StoreImageAllROIs(:,:,j))-StoreImageMean)./StoreImageMean)*100;
            display(j/size(StoreImageAllROIs,3)*100)
        end
        
        %average
        m1=mean(StoreImageAllROIs,1);
        m2=mean(m1,2);
        %average over the time
        MeanImFluoOverTime = squeeze(m2);
        clear m1 m2 StoreImageAllROIs
        
        %resample TimeMeanImage
        imFluoRes = resample(MeanImFluoOverTime,ResampleParRobotFluo,1);
        
        %%%to save data
        % Mean Fluo in ROI i       
        imFluoRes_filt    = imFluoRes;
        
        dataGCamp.fluoROI = imFluoRes_filt - median(imFluoRes_filt);
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     

        
        
        %%%%%%%%%%%%%%%
        %% SAVE DATA %%
        %%%%%%%%%%%%%%%
        
        %%%%%%
        % time
        dataGCamp.t = t_robot;
        % force status
        dataGCamp.status = status_robot;
        % force x
        dataGCamp.fx = fx_robot;
        if ~isempty(pos_robot)
            %positionhold
            dataGCamp.pos   = pos_robot;
            %speed
            dataGCamp.speed = speed_robot;
            %acceleration
            dataGCamp.acc   = acc_robot;
        end
        %%%%%%
        
        %%%%%%
        %%% Info Trial %%%
        IndexStatusMENO1 = find(status_robot==-1);
        %punto di start di ogni prova
        WholeTrials_Index = [IndexStatusMENO1(find(diff(IndexStatusMENO1)>1)); IndexStatusMENO1(end)];
        % number of completed trials
        NumTrials = length(WholeTrials_Index);
        dataGCamp.InfoTrial.NumTrials = NumTrials;
        
        TrialsVector = zeros(length(status_robot),1);
        WholeTrials_Index = [1; WholeTrials_Index];
        for nt=1:length(WholeTrials_Index)
            if nt<length(WholeTrials_Index)
                TrialsVector(WholeTrials_Index(nt): WholeTrials_Index(nt+1)-1) =  nt;
            else
                TrialsVector(WholeTrials_Index(nt): end) =  nt;
            end
        end
        % division of the array based on the trials (1-> trial 1; 2-> trial 2 etc.)
        dataGCamp.InfoTrial.TrialsVector  = TrialsVector;
        
        IndexStatus3 = find(status_robot == 3);
        %inizio dei trials (status3) nel vettore IndexTrials3
        Trials_In = [1; find(diff(IndexStatus3)>1)+1];
        %fine dei trials (status3) nel vettore IndexTrials3
        Trials_End = [Trials_In(2:end)-1; length(IndexStatus3)];
        
        % start and end points of the status 3 of the trials
        dataGCamp.InfoTrial.Trials_Start_End_St3  = [IndexStatus3(Trials_In(1:NumTrials)) IndexStatus3(Trials_End(1:NumTrials))];
        %%%%%%
        
        %%%%%%
        %%% Info %%%
        % name animal
        dataGCamp.Info.Name = Animal_Name_choice;
        % day
        dataGCamp.Info.Date = TrialDay(1:2);
        % treatment
        if strfind(Animal_Name_choice,'stroke')
            dataGCamp.Info.Treat = 1;
        elseif strfind(Animal_Name_choice,'control')
            dataGCamp.Info.Treat = 0;
        end
        % frequency of sampling
        dataGCamp.Info.Fs = fs_robot;
        % frequency of sampling of original Fluo
        dataGCamp.Info.Fs_FluoImages = fs_fluo;
        
        % interval Fluo and Force Mean
        dataGCamp.Info.IntervalToFluoForceMean_PointFluoFreq  = IntervalToFluoMean;
        dataGCamp.Info.IntervalToFluoForceMean_sec            = IntervalToFluoMean/dataGCamp.Info.Fs_FluoImages;
        dataGCamp.Info.IntervalToFluoForceMean_PoinForceFreq  = (IntervalToFluoMean/dataGCamp.Info.Fs_FluoImages)*dataGCamp.Info.Fs;
        
        
        % ROIs
        dataGCamp.Info.ROI = ROI_name;
        dataGCamp.Info.ROI_PosAndSize = ROI_PosAndSize';
        %%%%%%
        
        %%%%%%
        filename_GCamp = ['dataMouseGCamp_',Animal_Name,'_',TrialDay(1:2)];
        save([CurrAnDayFolder,'\',filename_GCamp],'dataGCamp')
        %%%%%%
        
        
        %%%%%%%%%%%%%%%
        %% PLOT DATA %%
        %%%%%%%%%%%%%%%
        figPl = figure;
        plot(dataGCamp.t,dataGCamp.status/10,'r')
        hold on
        plot(dataGCamp.t,dataGCamp.pos/10,'g')
        plot(dataGCamp.t,dataGCamp.fx,'b')
        plot(dataGCamp.t,dataGCamp.fluoROI/20,'g')
        legend([{'Status','Position','Force'},ROI_name])
        Day_Fig_filename = 'Fig_whole_sig';
        saveas(figPl,[CurrAnDayFolder,'\',Day_Fig_filename],'fig');
        pause(1.5)
        close(figPl)
        
        %%%%
        
    end
end

        
        
        

