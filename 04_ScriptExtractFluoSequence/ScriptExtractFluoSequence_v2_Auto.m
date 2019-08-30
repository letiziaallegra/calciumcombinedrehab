%
% script to extract fluo frames corresponding to force peaks
% the extraction of the frames is based on the previous detection of both
% force and fluo peaks (stored in dataGcamp..._Par in the field PeaksPar_Fx_Fluo)
%
%
% SALVA I LONG
%


clear all
close all
clc


CurrDir = cd;
AnimalMainDir = ['/Volumes/ALE6TB_DESK/DATA/Imaging/Emilia/Rehab/MATLAB_DATA_FOLDERS/'];
AnimalMainDir = '/Volumes/ALE6TB_DESK/DATA/Imaging/Emilia/toxin/MATLAB/';
AnimalMainDir = '/Users/alessandro/Desktop/toxin/MATLAB';
AnimalMainDir = '/Users/alessandro/Desktop/180424_RehabOptogen/MATLAB/';
AnimalMainDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB/';
AnimalMainDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_GCaMP27-29';
AnimalMainDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_OR';
AnimalMainDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_CTRL_190122';
AnimalMainDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_CTRL_190117';

%%%%%%%%%%%
% ListAnimalTogether = {      'GCaMPChR2_7_control', 'GCaMPChR2_17_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
%                             'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke','GCaMPChR2_25_stroke', 'GCaMPChR2_26_stroke'...
%                             'GCaMP16_stroke_BoNT','GCaMP17_stroke_BoNT','GCaMP18_stroke_BoNT',...
%                             'GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT',...
%                             'GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT','GCaMPChR2_16_stroke_BoNT'};
                        
ListAnimalTogether = { 'GCampChR2_TOX1','GCampChR2_TOX2'};      
ListAnimalTogether = {  'GCamp_24_sani', 'GCamp_25_sani', 'GCamp_26_sani', 'GCamp_22_onlyrob', 'GCamp_23_onlyrob', 'GCamp_24_onlyrob',...
    'GCamp_25_onlyrob', 'GCamp_26_onlyrob'}; 
ListAnimalTogether = {  'GCamp_26_sani' }; 
ListAnimalTogether = {  'Tox1_toxin','Tox2_toxin','Tox3_toxin','Tox4_toxin','Tox5_toxin' }; 
ListAnimalTogether = {  'or19_robot','or20_robot','or21_robot','or23_optostim','or24_optostim',...
    'or25_optostim','or26_optostim+robot','or27_robot','or28_optostim+robot', 'or29_sham', 'or30_sham'};
ListAnimalTogether = {'or27_robot'};
ListAnimalTogether = { 'GCaMP-ChR2-22_stroke', 'GCaMP20_Ctrl', 'GCaMP21_Ctrl' , 'GCaMP22_Ctrl', ...
    'GCaMP22_robot', 'GCaMP23_Ctrl', 'GCaMP23_robot', 'GCaMP24_robot', 'GCaMP25_robot', 'GCaMP26_robot', ...
    'GCaMP27_Ctrl', 'GCaMP28_Ctrl', 'GCaMP29_Ctrl'};
ListAnimalTogether = { 'GCaMP23_robot', 'GCaMP24_robot', 'GCaMP25_robot', 'GCaMP26_robot', ...
    'GCaMP27_Ctrl', 'GCaMP28_Ctrl', 'GCaMP29_Ctrl'};
ListAnimalTogether = { 'GCaMP-ChR2-22_stroke'};
ListAnimalTogether = { 'GCaMP27_wk02_Ctrl', 'GCaMP27_wk03_Ctrl','GCaMP27_wk04_Ctrl','GCaMP28_wk02_Ctrl'...
,'GCaMP28_wk03_Ctrl','GCaMP28_wk04_Ctrl','GCaMP29_wk02_Ctrl','GCaMP29_wk03_Ctrl','GCaMP29_wk04_Ctrl'};
ListAnimalTogether = {'OR9_wk01_sham', 'OR10_wk01_sham', 'OR13_wk01_optostim','OR14_wk01_optostim'...
'OR15_wk01_optostim+robot','OR16_wk01_optostim+robot','OR17_wk01_optostim+robot'};
ListAnimalTogether = {'GCaMP24_wk01_ctrl','GCaMP25_wk01_ctrl','GCaMP26_wk01_ctrl'};
ListAnimalTogether = {'GCaMP-ChR2-1_wk01_ctrl'};

%size of the cranial window (mm)
Size_CW = 4.4; %(mm)
% Size_CW = 5.25; %(mm)


%making the new folder structure make sure the following folder is the one
%used in the next script


for lat=1:length(ListAnimalTogether)
    
    SAVE_IMAGES_STACK = 1;
    
%     SAVE_IMAGES_STACK = 0;
    
    %%%%%%%%%%%
    Animal_name = ListAnimalTogether{lat};
    %%%%%%%%%%%
    
    %%%%  Folder
    UsbPort = 'M';
    AnimalDir     = ['C:\LENS\Data Leti'];
    AnimalDir = ['/Volumes/ALE6TB_DESK/DATA/Imaging/Emilia/Rehab/MATLAB_DATA_FOLDERS/'];
    AnimalDir = '/Volumes/ALE6TB_DESK/DATA/Imaging/Emilia/toxin/MATLAB/';
    AnimalDir = '/Users/alessandro/Desktop/toxin/MATLAB';
    AnimalDir = '/Users/alessandro/Desktop/180424_RehabOptogen/MATLAB/';
    AnimalDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB/';
    AnimalDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_GCaMP27-29';
    AnimalDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_OR';
    AnimalDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_CTRL_190122';
    AnimalDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_CTRL_190117';
    
    
    % AnimalDir     = [UsbPort,':\LENS\Animals Data\NoBregmaREF'];
    MainDir       = [AnimalDir,filesep,Animal_name];
    
    
    save_to_working_folder = 1;
    WORKING_FOLDER = ['/Volumes/ALE6TB_DESK/DATA/Imaging/Emilia/toxin/Sequence_Trials_DATA/'];
    WORKING_FOLDER = ['/Users/alessandro/Desktop/toxin/Sequence_Trials_DATA/'];
    WORKING_FOLDER = [fileparts(AnimalDir(1:end-1)),filesep,'Sequence_Trials_DATA/'];
     
    %%%%%%%% load Reference File %%%%%%%%%%%%%%%%%%%%%%%%
    RefDir      = ['C:\Users\asus\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\Script_Flip_Find_References\MAT_Rot_Trans'];
    RefDir = ['/Users/alessandro/Desktop/ELABORAZIONE DATA/Script_Flip_Find_References/MAT_Rot_Trans'];
    FileRefName = [Animal_name,'_Rot_Trans_Par.mat'];
    if exist([RefDir,filesep,FileRefName])
        % rot_transl
        load([RefDir,filesep,FileRefName]);
    else
        error([RefDir,filesep,FileRefName,' is not present in the folder']);
    end
    
    
    NumDaysFolder = dir(MainDir);
    
    for nd_i=3:length(NumDaysFolder)
        
        DayCurrDir = NumDaysFolder(nd_i,1).name;
        
        if isdir([MainDir,filesep,DayCurrDir])
            
            filename = ['dataMousegcamp_',Animal_name,'_',DayCurrDir(1:2),'_Par','.mat'];
            
            %%% Initial Check %%%
            if exist([MainDir,filesep,DayCurrDir,filesep,filename])
                
                %load file gCamp
                load([MainDir,filesep,DayCurrDir,filesep,filename])
                
                if isfield(dataGCamp,'PeaksPar_Fx_Fluo')
                    display('Check OK')
                else
                    error('PeaksPar_Fx_Fluo missing: perform the parameters extraction')
                end
            else
                error('load dataGCamp_Par')
            end
            %%%
            
            
            %%% Time and Frequencies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            t      = dataGCamp.t;
            Fs     = dataGCamp.Info.Fs;
            FsFluo = dataGCamp.Info.Fs_FluoImages; %[Hz]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%% Select Frame used as reference and ROI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            CentralFrame = [1 3];
            %%% ROI to use for peaks
            ROISelected  = 1; %[whole]
            %%% Number of Frame to Take around the CentralFrame
            FsFluo = dataGCamp.Info.Fs_FluoImages;
            BeforeCentralFrame = 0.32*FsFluo; %(0.32sec -> 8 frames)
            AfterCentralFrame  = 2*FsFluo; %(2sec)
            %%%
            downsamplingfactor = 1; %if FsFluo == 25
            % downsamplingfactor = 12; %if FsFluo == 100
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%% Info images store %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            Animal_Name = dataGCamp.Info.Name;
            TrialDay    = dataGCamp.Info.Date;
            %size ROIs
            ROI_PosAndSize  = dataGCamp.Info.ROI_PosAndSize;
            ROI_name        = dataGCamp.Info.ROI;
            %interval of frames to compute f0
            IntervalToFluoMean = dataGCamp.Info.IntervalToFluoForceMean_PointFluoFreq;
            %region of interest to take
            ROI_num  = length(ROI_name);
            %%%
            
            
            %%% data Images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isdir([MainDir,filesep,DayCurrDir,filesep,'MAT_trial'])
                folderTASK_FLUO = [MainDir,filesep,DayCurrDir,filesep,'MAT_trial'];
            else
                error([MainDir,filesep,DayCurrDir,filesep,'MAT_trial',' folder does not exist'])
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%% Images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            
            
            %%% Find correspondance between CentralFrame and Frames of the video %%%%%
            TimeTrigger = dataGCamp.PeaksPar_Fx_Fluo{CentralFrame(1),ROISelected}(:,CentralFrame(2)); %[100 Hz]
            StoreTrigger_Fluo = [];
            for tt=1:length(TimeTrigger)
                
                Cent     =  round((TimeTrigger(tt)*FsFluo));
                BefCent  =  Cent-BeforeCentralFrame;
                AftCent  =  Cent+AfterCentralFrame;
                
                StoreTrigger_Fluo = [StoreTrigger_Fluo; [BefCent Cent AftCent]];
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            index=0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%% find f0 on the ROI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            rw = ROI_PosAndSize(3,ROISelected)+1;
            cl = ROI_PosAndSize(4,ROISelected)+1;
            
            MatrixImageF0 = zeros(rw,cl,IntervalToFluoMean(2)-IntervalToFluoMean(1)+1);
            indexStore = 0;
            
            for fi=IntervalToFluoMean(1):IntervalToFluoMean(2)
                
                %index to report to images list
                indexImage = fi+2;
                %name
                nameImage = d(indexImage,1).name;
                
                %load image
                load([folderTASK_FLUO,filesep,nameImage]);
                Im = Im8_fv;
                Im_Original = Im;
                
                %filtering
                Im = medfilt2(Im);
                
                %
                indexStore = indexStore+1;
                MatrixImageF0(:,:,indexStore) = Im;
                
            end
            MeanMatrix_F0 = round(mean( double(MatrixImageF0),3)); %in gray tones
            MeanMatrix_F0(MeanMatrix_F0==0) = 1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%% find f0 on the ROI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            MatrixImageSequence = cell(size(StoreTrigger_Fluo,1),1);
            
            for stf=1:size(StoreTrigger_Fluo,1) %scroll the list of trigger
                
                Bef_buf = StoreTrigger_Fluo(stf,1);
                Cen_buf = StoreTrigger_Fluo(stf,2);
                Aft_buf = StoreTrigger_Fluo(stf,3);
                
                rw = ROI_PosAndSize(3,ROISelected)+1;
                cl = ROI_PosAndSize(4,ROISelected)+1;
                
                %List frames to take
                ListFramesToTake = [Bef_buf:downsamplingfactor:Aft_buf];
                
                MatrixImageSequence_buf = zeros(rw,cl,length(ListFramesToTake));
                
                indexStore = 0;
                
                for iL=1:length(ListFramesToTake)
                    
                    i = ListFramesToTake(iL);
                    
                    %index to report to images list
                    indexImage = i+2;
                    
                    %check number of frames (chosen interval goes over the limit)
                    if indexImage - length(d)>0
                        Im = ones(ROI_PosAndSize(3,ROISelected), ROI_PosAndSize(4,ROISelected))*NaN;
                    else
                        %name
                        nameImage = d(indexImage,1).name;
                        
                        %load image
                        load([folderTASK_FLUO,filesep,nameImage]);
                        Im = Im8_fv;
                        Im_Original = Im;
                        
                        %filtering
                        Im = medfilt2(Im);
                        
                        Im = double(Im);
                        
                        %delta_f/f0
                        Im_norm = ((Im-MeanMatrix_F0)./MeanMatrix_F0)*100;
                        
                        %CUT EDGES (it could be possible to have fake gray tones, put after the registration of the image)
                        minEdge = min(min(abs(Im_norm)));
                        Im_norm(:,[1:2])   = minEdge;
                        Im_norm(:,end)     = minEdge;
                        Im_norm(:,end-1)   = minEdge;
                        Im_norm([1 2],:)   = minEdge;
                        Im_norm(end,:)     = minEdge;
                        Im_norm(end-1,:)   = minEdge;
                        
                        
                    end
                    
                    %store
                    indexStore = indexStore+1;
                    MatrixImageSequence_buf(:,:,indexStore) = Im_norm;
                    
                    %%%%%%MAX AND MIN
                    if indexStore == 1
                        max_Im_norm_old = max(max(Im_norm));
                        min_Im_norm_old = min(min(Im_norm));
                    else
                        max_Im_norm = max(max(Im_norm));
                        min_Im_norm = min(min(Im_norm));
                        if max_Im_norm > max_Im_norm_old
                            max_Im_norm_old = max_Im_norm;
                        end
                        
                        if min_Im_norm < min_Im_norm_old
                            min_Im_norm_old = min_Im_norm;
                        end                            
                    end
                    %%%%%%%                   
                end
                
                MatrixImageSequence{stf,1} = MatrixImageSequence_buf;
                
            end
            
            %%% DELETE AREA AROUND BREGMA %%%%%%%%%%%%%%%%%%%%%%%%%
            for ii=1:size(MatrixImageSequence,1)
                MatrixImage=MatrixImageSequence{ii,1};
                for jj=1:size(MatrixImage,3)
                    i_day_actual = find(rot_transl(:,1) == str2num(DayCurrDir(1:2)));
                    dim=size(MatrixImage,2);
                    MatrixImage(dim-rot_transl(i_day_actual,5)-8:dim-rot_transl(i_day_actual,5)+8,:,jj)=0;
                end
                MatrixImageSequence{ii,1}=MatrixImage;
            end
            
            %%% SAVE MAT MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ImageSequence.Name                = dataGCamp.Info.Name;
            ImageSequence.Date                = dataGCamp.Info.Date;
            ImageSequence.MatrixImageSequence = MatrixImageSequence;
            ImageSequence.CentralFrameInfo    = CentralFrame;
            ImageSequence.ROISelected         = ROISelected;
            ImageSequence.CentFrame           = BeforeCentralFrame+1;
            ImageSequence.BeforeCentralFrame  = BeforeCentralFrame;
            ImageSequence.AfterCentralFrame   = AfterCentralFrame;
            ImageSequence.downsamplingfactor  = downsamplingfactor;
            
            Folder_Seq_Name = [MainDir,filesep,DayCurrDir,filesep,'SequenceLong'];
            if ~isdir(Folder_Seq_Name)
                
                %make Folder where saving
                mkdir(Folder_Seq_Name)
            end
            
            Filename = ['ImageSequence_',dataGCamp.Info.Name,'_',dataGCamp.Info.Date,'_CenFrame_',num2str(CentralFrame(1)),'_',num2str(CentralFrame(2)),'_ROI_',num2str(ROISelected)];
            save([Folder_Seq_Name,filesep,Filename],'ImageSequence','-v7.3');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if save_to_working_folder
                %additional copy structure for the next script
                if ~isdir(WORKING_FOLDER)
                    
                    %make Folder where saving
                    mkdir(WORKING_FOLDER)
                end
                new_animal_folder = [WORKING_FOLDER,filesep,Animal_Name];
                
                if ~isdir(new_animal_folder)
                    mkdir(new_animal_folder)
                end
                
                new_day_folder = [new_animal_folder,filesep,Animal_Name,'_',TrialDay,'_SequenceLong_TIF'];
                
                if ~isdir(new_day_folder)
                    mkdir(new_day_folder)
                end
                
                sequence_folder = [new_day_folder,filesep,'SequenceLong'];
                
                if ~isdir(sequence_folder)
                    mkdir(sequence_folder)
                end
                
                save([sequence_folder,filesep,Filename],'ImageSequence','-v7.3');
            end
            
            %%%%%%% save images %%%%%%%%%%%%
            if SAVE_IMAGES_STACK == 1
                
                Folder_Seq_NameTIF = [MainDir,filesep,DayCurrDir,filesep,[ImageSequence.Name ,'_',ImageSequence.Date,'_SequenceLong_TIF']];
                if ~isdir(Folder_Seq_NameTIF)
                    
                    %make Folder where saving
                    mkdir(Folder_Seq_NameTIF)
                end
                
                
                for mc=1:size(MatrixImageSequence,1)
                    
                    lenMIS = size(MatrixImageSequence{mc,1},3);
                    
                    indexImToSave = 10;
                    for mr=1:lenMIS
                        
                        indexImToSave    = indexImToSave+1;
                        ImToSave         = squeeze(ImageSequence.MatrixImageSequence{mc,1}(:,:,mr));
                        filenameImToSave = [ImageSequence.Name,'_',ImageSequence.Date,'_Seq_Image_',num2str(mc),'_',num2str(indexImToSave)];
                        %scaling
                        ImToSave_SCALED = (ImToSave-min_Im_norm_old) * (2^8-1-0)/(max_Im_norm_old-min_Im_norm_old);

                        imwrite(uint8( ImToSave_SCALED),[Folder_Seq_NameTIF,filesep,filenameImToSave,'.tif']);
                    
                    end
                   
                end
                clear ImToSave filenameImToSave ImToSave_SCALED
                FilenameText_ToSave = ['_Info_',ImageSequence.Name,'_',ImageSequence.Date,'_Seq_Image'];
                Text_ToSave         = [min_Im_norm_old max_Im_norm_old];
                save([Folder_Seq_NameTIF,filesep,FilenameText_ToSave,'.txt'],'Text_ToSave','-ASCII');
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            display(['End Process for', Filename])
            clear MatrixImageSequence CentralFrame ROISelected BeforeCentralFrame AfterCentralFrame
        else
            display([MainDir,filesep,DayCurrDir,' is not a directory'])
            
        end
    end
    
    close all
    clearvars -except lat ListAnimalTogether CurrDir
    
end
display('End Process')


