%
% script to extract fluo frames corresponding to force peaks
% the extraction of the frames is based on the previous detection of both
% force and fluo peaks (stored in dataGcamp..._Par in the field PeaksPar_Fx_Fluo)
% 

clear all
close all
clc


CurrDir = cd;

% ListAnimalTogether = {      'GCaMPChR2_17_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
%                             'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke','GCaMPChR2_25_stroke',...
%                             'GCaMP16_stroke_BoNT','GCaMP17_stroke_BoNT','GCaMP18_stroke_BoNT',...
%                             'GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT',...
%                             'GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT','GCaMPChR2_16_stroke_BoNT'};

ListAnimalTogether = {'GCaMPChR2_25_stroke_LateRehab'};

for lat=1:length(ListAnimalTogether)
%%%%%%%%%%%
Animal_name = ListAnimalTogether{lat};
%%%%%%%%%%%

%%%%  Folder
UsbPort = 'F';
AnimalDir     = [UsbPort,':\LENS\Animals Data'];
% AnimalDir     = [UsbPort,':\LENS\Animals Data\NoBregmaREF'];
MainDir       = [AnimalDir,'\',Animal_name];

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
    t         = dataGCamp.t;
    status    = dataGCamp.status;
    Fs        = dataGCamp.Info.Fs;
    FsFluo    = dataGCamp.Info.Fs_FluoImages; %[Hz]
    NumTrials = dataGCamp.InfoTrial.NumTrials;
    %%% ROI to use for peaks
    ROISelected  = 1; %[whole]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Select Frame used as reference and ROI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    status_REST_Selected = [5:7];
    status_REST = zeros(1,length(status));
    status_REST(status==status_REST_Selected(1) | status==status_REST_Selected(2)|status==status_REST_Selected(3)) = 1;
    d_status_REST = diff(status_REST);
       
    status_REST_St = find(d_status_REST==1)+1;
    status_REST_En = find(d_status_REST==-1);
    
    status_REST_St = status_REST_St(1:NumTrials);
    status_REST_En = status_REST_En(1:NumTrials);
    
    status_REST_St_sec = status_REST_St/Fs;
    status_REST_En_sec = status_REST_En/Fs;
    
    %minimum size of sequence
    MinSize = min(abs((status_REST_En_sec*FsFluo-status_REST_St_sec*FsFluo)));

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
    StoreTrigger_Fluo = [];
    for tt=1:NumTrials
        
        BefCent  =  round(status_REST_St_sec(tt)*FsFluo);
        AftCent  =  round(status_REST_En_sec(tt)*FsFluo);
        
        StoreTrigger_Fluo = [StoreTrigger_Fluo; [BefCent AftCent]];
        
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
        Aft_buf = StoreTrigger_Fluo(stf,2);
        
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
                                
%                 %CUT EDGES (it could be possible to have fake gray tones, put after the registration of the image)
%                 minEdge = min(min(abs(Im_norm)));
%                 Im_norm(:,[1:2])   = minEdge;
%                 Im_norm(:,end)     = minEdge;
%                 Im_norm(:,end-1)   = minEdge;
%                 Im_norm([1 2],:)   = minEdge;
%                 Im_norm(end,:)     = minEdge;
%                 Im_norm(end-1,:)   = minEdge;
                
                
            end
            
            %store
            indexStore = indexStore+1;
            MatrixImageSequence_buf(:,:,indexStore) = Im_norm;
            
            
        end
        
        MatrixImageSequence{stf,1} = MatrixImageSequence_buf;
        
    end
    
    
    %%% SAVE MAT MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ImageSequence.Name                 = dataGCamp.Info.Name;
    ImageSequence.Date                 = dataGCamp.Info.Date; 
    ImageSequence.MinSize              = floor(MinSize); 
    ImageSequence.MatrixImageSequence  = MatrixImageSequence;
    ImageSequence.status_REST_Selected = status_REST_Selected;
    ImageSequence.ROISelected          = ROISelected;
    ImageSequence.downsamplingfactor   = downsamplingfactor;
    
    Folder_Seq_Name = [MainDir,'\',DayCurrDir,'\','Sequence_REST'];
    if ~isdir(Folder_Seq_Name)
        
        %make Folder where saving
        mkdir(Folder_Seq_Name)
    end
    
    Filename = ['ImageSequence_',dataGCamp.Info.Name,'_',dataGCamp.Info.Date,'_Status_',num2str(status_REST_Selected),'_ROI_',num2str(ROISelected)];
    save([Folder_Seq_Name,'\',Filename],'ImageSequence','-v7.3');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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