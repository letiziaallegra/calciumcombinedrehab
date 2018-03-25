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
%%%%%%%%%%%
% ListAnimalTogether = {      'GCaMPChR2_7_control', 'GCaMPChR2_17_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
%                             'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke','GCaMPChR2_25_stroke', 'GCaMPChR2_26_stroke'...
%                             'GCaMP16_stroke_BoNT','GCaMP17_stroke_BoNT','GCaMP18_stroke_BoNT',...
%                             'GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT',...
%                             'GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT','GCaMPChR2_16_stroke_BoNT'};
                        
ListAnimalTogether = { 'GCampChR2_TOX1','GCampChR2_TOX2'};      

%size of the cranial window (mm)
%Size_CW = 4.4; %(mm)
Size_CW = 5.25; %(mm)

for lat=1:length(ListAnimalTogether)
    
    SAVE_IMAGES_STACK = 1;
%     SAVE_IMAGES_STACK = 0;
    
    %%%%%%%%%%%
    Animal_name = ListAnimalTogether{lat};
    %%%%%%%%%%%
    
    %%%%  Folder
    UsbPort = 'M';
    AnimalDir     = ['C:\LENS\Data Leti'];
    % AnimalDir     = [UsbPort,':\LENS\Animals Data\NoBregmaREF'];
    MainDir       = [AnimalDir,'\',Animal_name];
    
    %%%%%%%% load Reference File %%%%%%%%%%%%%%%%%%%%%%%%
    RefDir      = ['C:\Users\asus\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\Script_Flip_Find_References\MAT_Rot_Trans'];
    FileRefName = [Animal_name,'_Rot_Trans_Par.mat'];
    if exist([RefDir,'\',FileRefName])
        % rot_transl
        load([RefDir,'\',FileRefName]);
    else
        error([RefDir,'\',FileRefName,' is not present in the folder']);
    end
    
    
    NumDaysFolder = dir(MainDir);
    
    for nd_i=3:length(NumDaysFolder)
        
        DayCurrDir = NumDaysFolder(nd_i,1).name;
        
        if isdir([MainDir,'\',DayCurrDir])
            
            filename = ['dataMouseGCamp_',Animal_name,'_',DayCurrDir(1:2),'_Par','.mat'];
            
            %%% Initial Check %%%
            if exist([MainDir,'\',DayCurrDir,'\',filename])
                
                %load file gCamp
                load([MainDir,'\',DayCurrDir,'\',filename])
                
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
            if isdir([MainDir,'\',DayCurrDir,'\','MAT_trial'])
                folderTASK_FLUO = [MainDir,'\',DayCurrDir,'\','MAT_trial'];
            else
                error([MainDir,'\',DayCurrDir,'\','MAT_trial',' folder does not exist'])
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
                load([folderTASK_FLUO,'\',nameImage]);
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
                        load([folderTASK_FLUO,'\',nameImage]);
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
            
            Folder_Seq_Name = [MainDir,'\',DayCurrDir,'\','SequenceLong'];
            if ~isdir(Folder_Seq_Name)
                
                %make Folder where saving
                mkdir(Folder_Seq_Name)
            end
            
            Filename = ['ImageSequence_',dataGCamp.Info.Name,'_',dataGCamp.Info.Date,'_CenFrame_',num2str(CentralFrame(1)),'_',num2str(CentralFrame(2)),'_ROI_',num2str(ROISelected)];
            save([Folder_Seq_Name,'\',Filename],'ImageSequence','-v7.3');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%% save images %%%%%%%%%%%%
            if SAVE_IMAGES_STACK == 1
                
                Folder_Seq_NameTIF = [MainDir,'\',DayCurrDir,'\',[ImageSequence.Name ,'_',ImageSequence.Date,'_SequenceLong_TIF']];
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

                        imwrite(uint8( ImToSave_SCALED),[Folder_Seq_NameTIF,'\',filenameImToSave,'.tif']);
                    
                    end
                   
                end
                clear ImToSave filenameImToSave ImToSave_SCALED
                FilenameText_ToSave = ['_Info_',ImageSequence.Name,'_',ImageSequence.Date,'_Seq_Image'];
                Text_ToSave         = [min_Im_norm_old max_Im_norm_old];
                save([Folder_Seq_NameTIF,'\',FilenameText_ToSave,'.txt'],'Text_ToSave','-ASCII');
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            display(['End Process for', Filename])
            clear MatrixImageSequence CentralFrame ROISelected BeforeCentralFrame AfterCentralFrame
        else
            display([MainDir,'\',DayCurrDir,' is not a directory'])
            
        end
    end
    
    close all
    clearvars -except lat ListAnimalTogether CurrDir
    
end
display('End Process')