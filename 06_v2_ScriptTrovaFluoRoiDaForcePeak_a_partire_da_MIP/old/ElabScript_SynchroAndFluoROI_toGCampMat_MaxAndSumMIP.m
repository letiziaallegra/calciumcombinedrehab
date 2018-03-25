%
% script to syncrhonize fluo and force signals and
% define the Fluo ROIs signals
%
clear
close all
clc

CurrDir = cd;


%% Choice of the animal and trial day
UsbPortHD = 'I';
Animal_name = 'GCaMPChR2_7_control';

%%%%%%%%%%% dataGCamp (ROI arbitrarie) -> folder where saving data
MainDir = [UsbPortHD,':\LENS\_data_MAT_GCamp\'];


%%%%% LOAD data MIP %%%%%%%%%%%%%%%%%%%%
cd('_data_MIP_GCamp')
cd(Animal_name)
%RegionBound_Index
load([Animal_name,'_MIP_max_sum_RegionBound_Index']);

%rot_transl (MAX)
RotTranslFolder   = [UsbPortHD,':\LENS\Script_Flip_Find_References\MAT_Rot_Trans'];
RotTranslFilename = [Animal_name,'_Rot_Trans_Par'];
load( [RotTranslFolder,'\',RotTranslFilename];

load([Animal_name,'_max_Rot_Trans_Par']);
rot_transl_MaxSum{1,:}= rot_transl;
clear rot_transl
%rot_transl (MAX)
load([Animal_name,'_sum_Rot_Trans_Par']);
rot_transl_MaxSum{2,:} = rot_transl;
clear rot_transl
cd(CurrDir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ListDays = {'01','02','03','04','05'};
%size ROIs
ROI_name = {'MaxMIP','SumMIP'};
ROI_num  = length(ROI_name);
%%

for i_day = 1:length(ListDays)
    
    TrialDay_str     = ListDays{i_day};
    TrialDay_int = str2double(TrialDay_str);
    
    %load single dataGCamp matrix
    DirGCampMatrix      = [MainDir,'\',Animal_name,'\',TrialDay_str];
    FileNameGCampMatrix = ['dataMouseGCamp_',Animal_name,'_',TrialDay_str,'_Par'];
    
    %load matrix GCamp
    load([DirGCampMatrix,'\',FileNameGCampMatrix ]);
    
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
    UsbPort = UsbPortHD;
    %
    Animal_name = dataGCamp.Info.Name;
    TrialDay    = dataGCamp.Info.Date;
    %interval of frames to compute f0
    IntervalToFluoMean = dataGCamp.Info.IntervalToFluoForceMean_PointFluoFreq;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% data folder %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd ..
    cd('_load_animal_fun')
    [folderTASK_FLUO ~] = fun_loadGCaMPanimal(UsbPort,Animal_name,TrialDay);
    cd(CurrDir)
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

    for i=RealStart_Fluo :downsamplingfactor:   LenFLUO_Im
        
        index = index+1;
        indexImage = i+2;
        waitbar(indexImage/(LenFLUO_Im),wb)
        
        nameImage = d(indexImage,1).name;
        
        Im_Original = imread([folderTASK_FLUO,'\',nameImage]);
        
        rw        = size(Im_Original,1);
        cl        = size(Im_Original,2);
        
        
        for nR = 1:ROI_num
            
            i_day_actual = find(rot_transl_MaxSum{nR,1}(:,1) == TrialDay_int);
            
            %%%%%%%%% ROTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            degree = rot_transl_MaxSum{nR,1}(i_day_actual,2);
            if  degree~= 0
                Im_R = imrotate(Im_Original,degree,'crop');
            end
            
            Im_OR          = zeros(rw,cl);
            Im_OR_RotTrasl = zeros(rw,cl);
            
            %%%%%%%%% TRANSLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %left/right transl
            trans_lr  = rot_transl_MaxSum{nR,1}(i_day_actual,3);
            %translation Y along rows     
            trans_Y   = rot_transl_MaxSum{nR,1}(i_day_actual,4);
            %translation X along columns
            trans_X   = rot_transl_MaxSum{nR,1}(i_day_actual,5);
            
            %translation along rows           
            if trans_lr>0
                Im_OR(1:rw-trans_Y+1,:) = Im_R(trans_Y:end,:);
                trans_lr = -1;
            elseif trans_lr<0
                Im_OR(trans_Y:end,:) = Im_R(1:rw-trans_Y+1,:);
                trans_lr = +1;
            end
            
            %translation along columns
            Im_OR_RotTrasl(:,1:cl-trans_X+1) = Im_OR(:,trans_X:end);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
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
                    Im_f0 = imread([folderTASK_FLUO,'\',nameImage_f0]);
                    Im_Original_f0 = Im_f0;
                    
                    %filtering
                    Im_Original_f0 = medfilt2(Im_Original_f0);
                    
                    % Use the mask to select part of the image
                    Im_f0 = double(Im_Original_f0) .* MASK;
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
                H_Im_Round_Fig_filename = ['Fig_ROI_around_ROI_MIP_',ROI_name{nR},'_',Animal_name,'_',TrialDay_str];
                saveas(H_Im_Round_Fig ,[DirGCampMatrix,'\',H_Im_Round_Fig_filename],'fig');
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
                
        end
        
        
    end
    
    
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
    H_MCA_fig_filename = ['Fig_Sig_MIP_',Animal_name,'_',TrialDay_str];
    saveas(H_MCA_fig,[DirGCampMatrix,'\',H_MCA_fig_filename],'fig');
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
    dataGCamp.ROI_MIP.rot_transl_MaxSum  = rot_transl_MaxSum;
    % signals
    dataGCamp.ROI_MIP.ROI_Signal         = imFluoResMean;
    
    
    %%%%%%
    filename_GCamp = ['dataMouseGCamp_',dataGCamp.Info.Name,'_',dataGCamp.Info.Date,'_Par','_MIP' ];
    save([DirGCampMatrix,'\',filename_GCamp],'dataGCamp')
    %%%%%%
    
    display(['END PROCESS for: ',filename_GCamp]);
    clear dataGCamp MatrixImage_f0 MatrixImage_f0_ROI MeanMatrix_f0 Im Cx Cy circle_image circlemask circle_image_f0 d imFluoResMean imFluoResMean_Buf t xgrid ygrid Im_Original Im_Original_f0 Im_f0
end









