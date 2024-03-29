%
% script to syncrhonize fluo and force signals and
% define the Fluo ROIs signals
%
clear
close all
clc

CurrDir = cd;


%% Choice of the animal and trial day
UsbPortHD = 'M';
UserName  = getenv('username');
%%%%%%%%%%% Animal Dir %%%%%%%%%%%%%%%%%
%AnimalDir = [UsbPortHD,':\LENS\Animals Data STIM'];
AnimalDir = ['/Volumes/ALE6TB_DESK/DATA/Imaging/Emilia/Rehab/MATLAB_DATA_FOLDERS/'];
AnimalDir = '/Users/alessandro/Desktop/toxin/MATLAB/';
AnimalDir = '/Users/alessandro/Desktop/180424_RehabOptogen/MATLAB/';
AnimalDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB/';
AnimalDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_GCaMP27-29/';
AnimalDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_OR/';
AnimalDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_CTRL_190122/';
AnimalDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_CTRL_190117/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Animal_Name %%%%%%%%%%%%%%%%%%%
% ListAnimalTogether = { 'GCaMPChR2_7_control'};


% ListAnimalTogether = {  'GCaMPChR2_7_control','GCaMPChR2_17_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
%                         'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke','GCaMPChR2_25_stroke','GCaMPChR2_26_stroke',...
%                         'GCaMP16_stroke_BoNT','GCaMP18_stroke_BoNT',...
%                         'GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT',...
%                         'GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT','GCaMPChR2_16_stroke_BoNT'};

ListAnimalTogether = {  'GCampChR2_TOX1' , 'GCampChR2_TOX2', 'GCampChR2_TOX3', 'GCampChR2_TOX4', 'GCampChR2_TOX5'};
ListAnimalTogether = {  'GCamp_26_sani' };
ListAnimalTogether = {  'Tox1_toxin','Tox2_toxin','Tox3_toxin','Tox4_toxin','Tox5_toxin' };
ListAnimalTogether = {  'or19_robot','or20_robot','or21_robot','or23_optostim','or24_optostim',...
    'or25_optostim','or26_optostim+robot','or27_robot','or28_optostim+robot', 'or29_sham', 'or30_sham'};
ListAnimalTogether = { 'GCaMP-ChR2-22_stroke', 'GCaMP20_Ctrl', 'GCaMP21_Ctrl' , 'GCaMP22_Ctrl', ...
    'GCaMP22_robot', 'GCaMP23_Ctrl', 'GCaMP23_robot', 'GCaMP24_robot', 'GCaMP25_robot', 'GCaMP26_robot', ...
    'GCaMP27_Ctrl', 'GCaMP28_Ctrl', 'GCaMP29_Ctrl'};
ListAnimalTogether = {'GCaMP20_Ctrl', 'GCaMP21_Ctrl' , 'GCaMP22_Ctrl', ...
    'GCaMP22_robot', 'GCaMP23_Ctrl', 'GCaMP23_robot', 'GCaMP24_robot', 'GCaMP25_robot', 'GCaMP26_robot', ...
    'GCaMP27_Ctrl', 'GCaMP28_Ctrl', 'GCaMP29_Ctrl'};
ListAnimalTogether = { 'GCaMP27_wk02_Ctrl', 'GCaMP27_wk03_Ctrl','GCaMP27_wk04_Ctrl','GCaMP28_wk02_Ctrl'...
,'GCaMP28_wk03_Ctrl','GCaMP28_wk04_Ctrl','GCaMP29_wk02_Ctrl','GCaMP29_wk03_Ctrl','GCaMP29_wk04_Ctrl'};
ListAnimalTogether = {'OR9_wk01_sham', 'OR10_wk01_sham', 'OR13_wk01_optostim','OR14_wk01_optostim'...
'OR15_wk01_optostim+robot','OR16_wk01_optostim+robot','OR17_wk01_optostim+robot'};
ListAnimalTogether = {'GCaMP24_wk01_ctrl','GCaMP25_wk01_ctrl','GCaMP26_wk01_ctrl'};
ListAnimalTogether = {'GCaMP-ChR2-1_wk01_ctrl'};

wba = waitbar(0,'Processing');

for LATo = 1:length(ListAnimalTogether)
    
    
    
    clearvars -except wba CurrDir UsbPortHD UserName AnimalDir ListAnimalTogether LATo
    
    Animal_Name_choice                    = [ListAnimalTogether{LATo}];
    wba = waitbar(0,wba,['Whole Process. Subject: ',Animal_Name_choice]);
    myString = findall(wba,'String',['Subject: ',Animal_Name_choice]);
    set(myString,'Interpreter','none');
    %%%%%%%% TrialDay %%%%%%%%%%%%%%%%%%%%%%%
    % TrialDay_choice = '01';
    TrialDay_choice = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% LOAD data MIP and SIP %%%%%%%%%%%%
    %ROI_MIP_Folder   = ['C:\Users\',UserName,'\Google Drive_SL\Piattaforma Stefano\ELABORAZIONE DATA\05_Maximum_Intensity_Projection_SeqSelez\_Store_MIP_SIP_RegionBound\',Animal_Name_choice];
    ROI_MIP_Folder   = ['/Users/alessandro/Desktop/ELABORAZIONE DATA/05_Maximum_Intensity_Projection_SeqSelez/_Store_MIP_SIP_RegionBound/',Animal_Name_choice];
    AnimalWhereSaveGCamp = ['/Users/alessandro/Desktop/ELABORAZIONE DATA/_data_MAT_GCamp_Store'];
    if ~isempty(strfind(Animal_Name_choice,'BoNT')) | strfind(Animal_Name_choice,'roboti')
        %load 4 files -> one a week
        % list of weeks !!!! IT IS IMPORTANT THAT THE WEEKDAY NUMBER IS
        % MATCHES THE REAL NUMBER, I.E. DAY #6 => FIRST DAY OF SECOND WEEK
        % Str_Week = {'Week_1_'; 'Week_2_'; 'Week_3_'; 'Week_4_'};
        % FOR ONE WEEK IT DOES NOT MATTER
        WeekList = {'Week_1_'; 'Week_2_'; 'Week_3_'; 'Week_4_'};
        weekSpace = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15; 16 17 18 19 20];
        WeekList = {'Week_1_'; 'Week_2_'; 'Week_3_'};
        weekSpace = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15];
    elseif strfind(Animal_Name_choice,'stroke_STIM_1w')
        WeekList = {''};
        weekSpace = [1:3];
    elseif strfind(Animal_Name_choice,'stroke_STIM')
        WeekList = {''};
        weekSpace = [16:20];
    else
        WeekList = {''};
        weekSpace = [1 2 3 4 5];
    end
    
    lenWeekList       = length(WeekList);
    RegionBound_Store = cell(lenWeekList,1);
    
    for i_wkl=1:lenWeekList
        
        ROI_MIP_Filename = [Animal_Name_choice,'_MIP_SIP_',WeekList{i_wkl},'RegionBound_Index'];
        %load RegionBound_Index
        load( [ROI_MIP_Folder,filesep,ROI_MIP_Filename] );
        %store RegionBound_Index for every week
        RegionBound_Store{i_wkl,1} =  RegionBound_Index;
    end
    %
    %ROIs (the same)
    ROI_name = RegionBound_Index(:,2)';
    ROI_num  = length(ROI_name);
    clear RegionBound_Index
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% rot_transl %%%%%%%%%%%%%%%%%%%%%%%
    %RotTranslFolder   = ['C:\Users\',UserName,'\Google Drive_SL\Piattaforma Stefano\ELABORAZIONE DATA\Script_Flip_Find_References\MAT_Rot_Trans'];
    RotTranslFolder   = ['C:\Users\asus\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\Script_Flip_Find_References\MAT_Rot_Trans'];
    RotTranslFolder = ['/Users/alessandro/Desktop/ELABORAZIONE DATA/Script_Flip_Find_References/MAT_Rot_Trans'];
    RotTranslFilename = [Animal_Name_choice,'_Rot_Trans_Par'];
    %load rot_transl
    load( [RotTranslFolder,filesep,RotTranslFilename] );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% choice of the animals to elaborate     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ListFolderAnimalDir = dir(AnimalDir);
    indexAn = 1;
    for lfcd=3:length(ListFolderAnimalDir)
        AnimalName_buf = ListFolderAnimalDir(lfcd,1).name;
        
        if length(AnimalName_buf)>5
            
            if strcmp(AnimalName_buf(1:2),'OR') & isdir([AnimalDir,filesep,AnimalName_buf])
                %all of the animals in folder
                AnimalName_IndexList(indexAn,1) = lfcd;
                indexAn = indexAn+1;
                
            end
            
            if strcmp(AnimalName_buf,Animal_Name_choice) & isdir([AnimalDir,filesep,AnimalName_buf])
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
    
    for anAn_i=1:length(AnimalName_IndexList)  %for animals
        
        Animal_Index = AnimalName_IndexList(anAn_i);
        Animal_Name  = ListFolderAnimalDir(Animal_Index,1).name;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% choice of the days to elaborate        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        AnimalCurrDir = dir([AnimalDir,filesep,Animal_Name]);
        indexDay = 1;
        
        for acd=3:length(AnimalCurrDir)
            AnimalDay_buf = AnimalCurrDir(acd,1).name;
            
            if isdir([AnimalDir,filesep,Animal_Name,filesep,AnimalDay_buf])
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
        for anD_i=1:length(Days_IndexList) %for scroll days
            
            TrialDay_Index = Days_IndexList(anD_i);
            TrialDay       = AnimalCurrDir(TrialDay_Index,1).name;
            
            CurrAnDayFolder = [AnimalDir,filesep,Animal_Name,filesep,TrialDay];
            CurrAnDayFolder_List = dir(CurrAnDayFolder);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%% find image folder and force file       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for cadf=3:length(CurrAnDayFolder_List)
                
                %additional check on name length
                if length((CurrAnDayFolder_List(cadf,1).name)) >2
                    if strcmp(CurrAnDayFolder_List(cadf,1).name(1:3),'MAT')
                        %data images folder
                        folderTASK_FLUO = [CurrAnDayFolder,filesep,CurrAnDayFolder_List(cadf,1).name];
                    elseif length(CurrAnDayFolder_List(cadf,1).name)>7
                        if strcmp(CurrAnDayFolder_List(cadf,1).name(end-7:end-4),'_Par')
                            %dataGCmap filename
                            folderTASK_dataGCampFilename = [CurrAnDayFolder,filesep,CurrAnDayFolder_List(cadf,1).name];
                        end
                    end
                end
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %load GCamp
            load(folderTASK_dataGCampFilename);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%% Initial Check %%%
            if exist('dataGCamp')
                if isfield(dataGCamp,'PeaksPar_Fx_Fluo')
                    display('Check OK')
                else
                    error('PeaksPar_Fx_Fluo missing: perform the parameters extraction')
                end
            else
                error('load dataGCamp_Par')
            end
            %%%
            
            
            %%% Info images store %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Animal_name = dataGCamp.Info.Name;
            TrialDay    = dataGCamp.Info.Date;
            %interval of frames to compute f0
            IntervalToFluoMean = dataGCamp.Info.IntervalToFluoForceMean_PointFluoFreq;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Time and Frequencies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            t      = dataGCamp.t;
            Fs     = dataGCamp.Info.Fs;
            FsFluo = dataGCamp.Info.Fs_FluoImages; %[Hz]
            ResampleParRobotFluo = Fs/FsFluo;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%% data Images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            d = dir(folderTASK_FLUO);
            
            len = length(d);
            NumImages = len-2;
            
            % - starting point
            RealStart_Fluo     = 1;
            
            %resampling fsample Robot Fluo
            ResampleParRobotFluo = round(Fs/FsFluo);
            %duration (it should be already ok because of the initial check done in 02_ScriptTrovaFluoRoiDaForcePeak)
            RealDurTask_Fluo = floor(length(t)/ResampleParRobotFluo);
            
            % - ending point
            LenFLUO_Im = RealStart_Fluo+RealDurTask_Fluo-1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            downsamplingfactor=1;
            
            index=0;
            wb_string = [Animal_Name_choice,' | ', num2str(TrialDay), '- Images Loading, Please wait...'];
            wb = waitbar(0, wb_string);
            myString = findall(wb, 'String', wb_string);
            set(myString,'Interpreter','none');
            wbpos = get(wba, 'Position');
            wbpos(2) = wbpos(2) + wbpos(4);
            set(wba, 'Position', wbpos)
            %%load image stack
            %% N.B. size of StoreImage depends on num of ROIs
            LenFLUO_Im = RealStart_Fluo+RealDurTask_Fluo-1;
            ROI_Centroid_MeanSignal = [];
            
            
            for i=RealStart_Fluo :downsamplingfactor:   LenFLUO_Im  %for images
                
                index = index+1;
                indexImage = i+2;
                
                
                nameImage = d(indexImage,1).name;
                
                %load .mat
                load([folderTASK_FLUO,filesep,nameImage]);
                Im = Im8_fv;
                Im_Original = Im;
                
                rw        = size(Im_Original,1);
                cl        = size(Im_Original,2);
                
                for nR = 1:ROI_num %for nR ROI
                    
                    i_day_actual = find(rot_transl(:,1) ==   str2num(TrialDay(1:2)) );
                    
                    %%%%%%%%% ROTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    degree = rot_transl(i_day_actual,2);
                    if  degree~= 0
                        Im_R = imrotate(Im_Original,degree,'crop');
                    else
                        Im_R = Im_Original;
                    end
                    
                    
                    Im_OR          = zeros(rw,cl);
                    Im_OR_RotTrasl = zeros(rw,cl);
                    
                    %%%%%%%%% TRANSLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %left/right transl
                    trans_lr  = rot_transl(i_day_actual,3);
                    %translation Y along rows
                    trans_Y   = rot_transl(i_day_actual,4);
                    %translation X along columns
                    trans_X   = rot_transl(i_day_actual,5);
                    
                    %translation along rows
                    if trans_lr<0
                        if round(trans_Y) == 0
                            trans_Y = 1;
                        end
                        Im_OR(1:round(rw-trans_Y+1),:) = Im_R(round(trans_Y):end,:);
                        trans_lr = -1;
                    elseif trans_lr>0
                        Im_OR(round(trans_Y):end,:) = Im_R(1:round(rw-trans_Y+1),:);
                        trans_lr = +1;
                    end
                    
                    %translation along columns
                    Im_OR_RotTrasl(:,1:cl-trans_X+1) = Im_OR(:,trans_X:end);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %Create Mask
                    if index==1
                        
                        %find the right week
                        [week_ok_index week_ok_index_r] =  find(weekSpace == str2num(TrialDay(1:2)));
                        
                        %Region Bounds
                        RegionBound_Index = RegionBound_Store{week_ok_index,1};
                        
                        MASK = zeros(rw,cl);
                        MASK(sub2ind([rw cl],RegionBound_Index{nR,1}(:,1),RegionBound_Index{nR,1}(:,2))) = 1;
                        MASK = imfill(MASK);
                        
                        %%% Trova f0 di questa area  %%%
                        MatrixImage_f0 = zeros(rw,cl,IntervalToFluoMean(2)-IntervalToFluoMean(1)+1);
                        indexStore_f0 = 0;
                        
                        for fi=IntervalToFluoMean(1):IntervalToFluoMean(2)
                            
                            %index to report to images list
                            indexImage_f0 = fi+2;
                            %name
                            nameImage_f0 = d(indexImage_f0,1).name;
                            
                            %load image
                            load([folderTASK_FLUO,filesep,nameImage_f0]);
                            Im_f0 = Im8_fv;
                            Im_Original_f0 = Im_f0;
                            
                            %filtering
                            Im_Original_f0 = medfilt2(Im_Original_f0);
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%% ROTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            degree = rot_transl(i_day_actual,2);
                            if  degree~= 0
                                Im_R_f0 = imrotate(Im_Original_f0,degree,'crop');
                            else
                                Im_R_f0 = Im_Original_f0;
                            end
                            
                            
                            Im_OR_f0          = zeros(rw,cl);
                            Im_OR_RotTrasl_f0 = zeros(rw,cl);
                            
                            %%%%%%%%% TRANSLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %left/right transl
                            trans_lr  = rot_transl(i_day_actual,3);
                            %translation Y along rows
                            trans_Y   = rot_transl(i_day_actual,4);
                            %translation X along columns
                            trans_X   = rot_transl(i_day_actual,5);
                            
                            %translation along rows
                            trans_Y = round(trans_Y);
                            if trans_lr<0
                                if round(trans_Y) == 0
                                    trans_Y = 1;
                                end
                                Im_OR_f0(1:rw-trans_Y+1,:) = Im_R_f0(trans_Y:end,:);
                                trans_lr = -1;
                            elseif trans_lr>0
                                Im_OR_f0(trans_Y:end,:) = Im_R_f0(1:rw-trans_Y+1,:);
                                trans_lr = +1;
                            end
                            
                            %translation along columns
                            Im_OR_RotTrasl_f0(:,1:cl-trans_X+1) = Im_OR_f0(:,trans_X:end);
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            % Use the mask to select part of the image
                            Im_f0 = double(Im_OR_RotTrasl_f0) .* MASK;
                            Im_f0(Im_f0==0) = NaN;
                            
                            %
                            indexStore_f0 = indexStore_f0+1;
                            MatrixImage_f0(:,:,indexStore_f0) = Im_f0;
                            
                        end
                        display('F0 extracted')
                        MeanMatrix_f0 = nanmean( double(MatrixImage_f0),3); %in gray tones
                        MeanMatrix_f0(MeanMatrix_f0==0) = 1;
                        
                        H_Im_Round_Fig = figure('Name',['ROI_around_',ROI_name{nR}]);
                        imshow(MeanMatrix_f0, []);
                        H_Im_Round_Fig_filename = ['Fig_ROI_around_',ROI_name{nR},'_',Animal_name,'_',TrialDay];
                        saveas(H_Im_Round_Fig ,[CurrAnDayFolder,filesep,H_Im_Round_Fig_filename],'fig');
                        close(H_Im_Round_Fig)
                        
                        %storage delle f0 delle ROI
                        MatrixImage_f0_ROI(:,:,nR) = MeanMatrix_f0;
                        
                        %storage delle f0 delle ROI
                        MatrixImage_MASK(:,:,nR) = MASK;
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    end
                    
                    waitbar(indexImage/(LenFLUO_Im),wb)
                    % Use the mask to select part of the image
                    Im_ROI = double(Im_OR_RotTrasl) .* MatrixImage_MASK(:,:,nR);
                    %Metto tutti gli 0 pari a NaN
                    Im_ROI(Im_ROI==0) = NaN;
                    
                    %fai delta_f/f0
                    Im = (Im_ROI-squeeze(MatrixImage_f0_ROI(:,:,nR)))./squeeze(MatrixImage_f0_ROI(:,:,nR))*100;
                    
                    %compute average of the image (just the circle bc the outside the circle is all NaN
                    ROI_MIP(i,nR)  = squeeze(nanmean(nanmean(Im,2),1));
                    
                    
                end %end for nR ROI
                
                
            end %end for images
            
            
            imFluoResMean = [];
            for nR = 1:length(ROI_name)
                
                imFluoResMean_Buf = resample(ROI_MIP(:,nR),ResampleParRobotFluo,1);
                imFluoResMean_Buf = imFluoResMean_Buf - median(imFluoResMean_Buf);
                
                imFluoResMean     = [imFluoResMean imFluoResMean_Buf];
            end
            
            
            
            
            %%%%%%%%%%%%%%%
            %% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%
            
            % image f0 ROI
            dataGCamp.ROI_MIP_SIP.ImageROI_f0        = MatrixImage_f0;
            % RegionBoundIndex
            dataGCamp.ROI_MIP_SIP.RegionBound_Index  = RegionBound_Index;
            % rot_transl_MaxSum
            dataGCamp.ROI_MIP_SIP.rot_transl         = rot_transl;
            % signals
            dataGCamp.ROI_MIP_SIP.ROI_Signal         = imFluoResMean;
            
            
            
            %%%%%%%%%%%%%%%
            %% PLOT DATA %%
            %%%%%%%%%%%%%%%
            H_MCA_fig = figure('Name','Fig_Sig_MIP_SIP_Area');
            hold on
            plot(dataGCamp.t,dataGCamp.status/10,'r');
            plot(dataGCamp.t,dataGCamp.fx,'b');
            plot(dataGCamp.t,dataGCamp.pos/10,'g');
            plot(dataGCamp.t,dataGCamp.speed/100,'m');
            plot(dataGCamp.t,dataGCamp.ROI_MIP_SIP.ROI_Signal(:,1)/10,'k');
            plot(dataGCamp.t,dataGCamp.ROI_MIP_SIP.ROI_Signal(:,2)/10,'g');
            plot(dataGCamp.t,dataGCamp.ROI_MIP_SIP.ROI_Signal(:,3)/10,'c');
            legend({'Status','Force','Pos/10','Speed/100','MIP/10','rostral SIP/10','caudal SIP/10'})
            H_MCA_fig_filename = ['Fig_Sig_MIP_SIP_',Animal_name,'_',TrialDay];
            saveas(H_MCA_fig,[CurrAnDayFolder,filesep,H_MCA_fig_filename],'fig');
            close(H_MCA_fig)
            %%%%
            
            %%%%%%
            filename_GCamp = ['dataMouseGCamp_',dataGCamp.Info.Name,'_',dataGCamp.Info.Date,'_Par','_MIPSIP' ];
            save([CurrAnDayFolder,filesep,filename_GCamp],'dataGCamp')
            dst_folder = [AnimalWhereSaveGCamp,filesep,Animal_Name];
            if ~isdir(dst_folder)
                mkdir(dst_folder)
            end
            save([dst_folder,filesep,filename_GCamp],'dataGCamp')
            %%%%%%
            
            display(['END PROCESS for: ',filename_GCamp]);
            clear dataGCamp MatrixImage_f0 MatrixImage_f0_ROI MeanMatrix_f0 Im Cx Cy circle_image circlemask circle_image_f0 d imFluoResMean imFluoResMean_Buf t xgrid ygrid Im_Original Im_Original_f0 Im_f0 ROI_MIP RegionBound_Index
            
        close(wb)
        end %end for scroll days
        
    waitbar(LATo/ length(ListAnimalTogether), wba, '')    
    end %end for animals
    
end %end LATo
close(wba)






