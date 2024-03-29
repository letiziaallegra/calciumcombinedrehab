%
% script to syncrhonize fluo and force signals and
% define the Fluo ROIs signals
%
clear
close all
clc

CurrDir = cd;


%% Choice of the animal and trial day
UsbPortHD = 'F';
UserName  = getenv('username');
%%%%%%%%%%% Animal Dir %%%%%%%%%%%%%%%%%
AnimalDir = [UsbPortHD,':\LENS\Animals Data\NoBregmaREF'];
% AnimalDir = ['L',':\LENS\Animals Data'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Animal_Name %%%%%%%%%%%%%%%%%%%
% ListAnimalTogether = { 'GCaMP3_control'};


ListAnimalTogether = {  ...
                        'GCaMP9_stroke','GCaMP10_stroke', 'GCaMP11_stroke', 'GCaMP14_stroke', 'GCaMP15_stroke',...
                        'GCaMPChR2_1_control'};



for LATo = 1:length(ListAnimalTogether)
    
    clearvars -except CurrDir UsbPortHD UserName AnimalDir ListAnimalTogether LATo
    
    Animal_Name_choice                    = [ListAnimalTogether{LATo}];
    
    %%%%%%%% TrialDay %%%%%%%%%%%%%%%%%%%%%%%
    % TrialDay_choice = '01';
    TrialDay_choice = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% LOAD data MIP %%%%%%%%%%%%%%%%%%%%
    ROI_MIP_Folder   = ['C:\Users\',UserName,'\Google Drive\Piattaforma Stefano\LENS_SSSA\ELABORAZIONE DATA\05_Maximum_Intensity_Projection\_Store_MIP_RegionBound\NoBregmaREF\',Animal_Name_choice];
    ROI_MIP_Filename = [Animal_Name_choice,'_MIP_RegionBound_Index'];
    %load RegionBound_Index
    load( [ROI_MIP_Folder,'\',ROI_MIP_Filename] );
    %
    %ROIs
    ROI_name = RegionBound_Index(:,2)';
    ROI_num  = length(ROI_name);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% rot_transl %%%%%%%%%%%%%%%%%%%%%%%
    RotTranslFolder   = ['C:\Users\',UserName,'\Google Drive\Piattaforma Stefano\LENS_SSSA\ELABORAZIONE DATA\Script_Flip_Find_References\MAT_Rot_Trans_NoBregmaREF'];
    RotTranslFilename = [Animal_Name_choice,'_Rot_Trans_Par'];
    %load rot_transl
    load( [RotTranslFolder,'\',RotTranslFilename] );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
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
    
    for anAn_i=1:length(AnimalName_IndexList)  %for animals
        
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
        for anD_i=1:length(Days_IndexList) %for scroll days
            
            TrialDay_Index = Days_IndexList(anD_i);
            TrialDay       = AnimalCurrDir(TrialDay_Index,1).name;
            
            CurrAnDayFolder = [AnimalDir,'\',Animal_Name,'\',TrialDay];
            CurrAnDayFolder_List = dir(CurrAnDayFolder);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%% find image folder and force file       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for cadf=3:length(CurrAnDayFolder_List)
                
                if strcmp(CurrAnDayFolder_List(cadf,1).name(1:3),'MAT')
                    %data images folder
                    folderTASK_FLUO = [CurrAnDayFolder,'\',CurrAnDayFolder_List(cadf,1).name];
                elseif length(CurrAnDayFolder_List(cadf,1).name)>7
                    if strcmp(CurrAnDayFolder_List(cadf,1).name(end-7:end-4),'_Par')
                        %dataGCmap filename
                        folderTASK_dataGCampFilename = [CurrAnDayFolder,'\',CurrAnDayFolder_List(cadf,1).name];
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
            wb = waitbar(0,'Images Loading, Please wait...');
            
            %%load image stack
            %% N.B. size of StoreImage depends on num of ROIs
            LenFLUO_Im = RealStart_Fluo+RealDurTask_Fluo-1;
            ROI_Centroid_MeanSignal = [];
            
            
            for i=RealStart_Fluo :downsamplingfactor:   LenFLUO_Im  %for images
                
                index = index+1;
                indexImage = i+2;
                waitbar(indexImage/(LenFLUO_Im),wb)
                
                nameImage = d(indexImage,1).name;
                
                %load .mat
                load([folderTASK_FLUO,'\',nameImage]);
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
                    Im_OR_2        = zeros(rw,cl);
%                     Im_OR_RotTrasl = zeros(rw,cl);
                    
                    %%%%%%%%% TRANSLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %left/right transl
                    trans_lr  = rot_transl(i_day_actual,3);
                    %up/down transl
                    trans_ud  = rot_transl(i_day_actual,6);
                    %translation Y along rows
                    trans_Y   = rot_transl(i_day_actual,4);
                    %translation X along columns
                    trans_X   = rot_transl(i_day_actual,5);
                    
                    %translation along rows
                    if trans_lr<0
                        Im_OR(1:rw-trans_Y+1,:) = Im_R(trans_Y:end,:);
                        trans_lr = -1;
                    elseif trans_lr>0
                        Im_OR(trans_Y:end,:) = Im_R(1:rw-trans_Y+1,:);
                        trans_lr = +1;
                    elseif trans_lr == 0
                        Im_OR = Im_R;
                        trans_lr = 0;
                    end
                    
                    %translation along colums
                    if trans_ud<0
                        Im_OR_2(:,1:cl-trans_X+1) = Im_OR(:,trans_X:end);
                        trans_ud = -1;
                    elseif trans_ud>0
                        Im_OR_2(:,trans_X:end) = Im_OR(:,1:cl-trans_X+1);
                        trans_ud = +1;
                    elseif trans_ud == 0
                        Im_OR_2 = Im_OR;
                        trans_ud = 0;
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    %final roto-translated image
                    Im_OR_RotTrasl = Im_OR_2;
                    
                    %Create Mask
                    if index==1
                        
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
                            load([folderTASK_FLUO,'\',nameImage_f0]);
                            Im_f0 = Im8_fv;
                            Im_Original_f0 = Im_f0;
                            
                            %filtering
                            Im_Original_f0 = medfilt2(Im_Original_f0);
                            
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%% ROTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            degree = rot_transl(i_day_actual,2);
                            if  degree~= 0
                                Im_R = imrotate(Im_Original_f0,degree,'crop');
                            else
                                Im_R = Im_Original_f0;
                            end
                            
                            
                            Im_OR_f0          = zeros(rw,cl);
                            Im_OR_2_f0        = zeros(rw,cl);
                            %                     Im_OR_RotTrasl = zeros(rw,cl);
                            
                            %%%%%%%%% TRANSLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %left/right transl
                            trans_lr  = rot_transl(i_day_actual,3);
                            %up/down transl
                            trans_ud  = rot_transl(i_day_actual,6);
                            %translation Y along rows
                            trans_Y   = rot_transl(i_day_actual,4);
                            %translation X along columns
                            trans_X   = rot_transl(i_day_actual,5);
                            
                            %translation along rows
                            if trans_lr<0
                                Im_OR_f0(1:rw-trans_Y+1,:) = Im_R_f0(trans_Y:end,:);
                                trans_lr = -1;
                            elseif trans_lr>0
                                Im_OR_f0(trans_Y:end,:) = Im_R_f0(1:rw-trans_Y+1,:);
                                trans_lr = +1;
                            elseif trans_lr == 0
                                Im_OR_f0 = Im_R_f0;
                                trans_lr = 0;
                            end
                            
                            %translation along colums
                            if trans_ud<0
                                Im_OR_2_f0(:,1:cl-trans_X+1) = Im_OR_f0(:,trans_X:end);
                                trans_ud = -1;
                            elseif trans_ud>0
                                Im_OR_2_f0(:,trans_X:end) = Im_OR_f0(:,1:cl-trans_X+1);
                                trans_ud = +1;
                            elseif trans_ud == 0
                                Im_OR_2_f0 = Im_OR_f0;
                                trans_ud = 0;
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            % Use the mask to select part of the image
                            Im_f0 = double(Im_OR_2_f0) .* MASK;
                            Im_f0(Im_f0==0) = NaN;
                            
                            %
                            indexStore_f0 = indexStore_f0+1;
                            MatrixImage_f0(:,:,indexStore_f0) = Im_f0;
                            
                        end
                        display('F0 extracted')
                        MeanMatrix_f0 = nanmean( double(MatrixImage_f0),3); %in gray tones
                        MeanMatrix_f0(MeanMatrix_f0==0) = 1;
                        
                        H_Im_Round_Fig = figure('Name',['ROI_around_ROI_MIP',ROI_name{nR}]);
                        imshow(MeanMatrix_f0, []);
                        H_Im_Round_Fig_filename = ['Fig_ROI_around_ROI_MIP_',ROI_name{nR},'_',Animal_name,'_',TrialDay];
                        saveas(H_Im_Round_Fig ,[CurrAnDayFolder,'\',H_Im_Round_Fig_filename],'fig');
                        close(H_Im_Round_Fig)
                        
                        %storage delle f0 delle ROI
                        MatrixImage_f0_ROI(:,:,nR) = MeanMatrix_f0;
                        
                        %storage delle f0 delle ROI
                        MatrixImage_MASK(:,:,nR) = MASK;
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    end
                    
                    
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
            %% PLOT DATA %%
            %%%%%%%%%%%%%%%
            H_MCA_fig = figure('Name','Fig_Sig_MIP_Area');
            plot(t,dataGCamp.status/10,'r')
            hold on
            plot(t,dataGCamp.fx,'b')
            plot(t,dataGCamp.speed/100,'m')
            plot(t,imFluoResMean/100)
            legend([{'Status','Force','Speed'},ROI_name])
            %saving figure
            H_MCA_fig_filename = ['Fig_Sig_MIP_',Animal_name,'_',TrialDay];
            saveas(H_MCA_fig,[CurrAnDayFolder,'\',H_MCA_fig_filename],'fig');
            close(H_MCA_fig)
            %%%%
            
            
            %%%%%%%%%%%%%%%
            %% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%
            
            % image f0 ROI
            dataGCamp.ROI_MIP.ImageROI_f0        = MatrixImage_f0;
            % RegionBoundIndex
            dataGCamp.ROI_MIP.RegionBound_Index  = RegionBound_Index;
            % rot_transl_MaxSum
            dataGCamp.ROI_MIP.rot_transl         = rot_transl;
            % signals
            dataGCamp.ROI_MIP.ROI_Signal         = imFluoResMean;
            
            
            %%%%%%
            filename_GCamp = ['dataMouseGCamp_',dataGCamp.Info.Name,'_',dataGCamp.Info.Date,'_Par','_MIP' ];
            save([CurrAnDayFolder,'\',filename_GCamp],'dataGCamp')
            %%%%%%
            
            display(['END PROCESS for: ',filename_GCamp]);
            clear dataGCamp MatrixImage_f0 MatrixImage_f0_ROI MeanMatrix_f0 Im Cx Cy circle_image circlemask circle_image_f0 d imFluoResMean imFluoResMean_Buf t xgrid ygrid Im_Original Im_Original_f0 Im_f0 ROI_MIP
            
            
        end %end for scroll days
        
        
    end %end for animals
    
end %end LATo







