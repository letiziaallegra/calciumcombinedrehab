%
% script to extract fluo frames corresponding to force peaks
% the extraction of the frames is based on the previous detection of both
% force and fluo peaks (stored in dataGcamp..._Par in the field PeaksPar_Fx_Fluo)
% 

% close all
clc
CurrDir = cd;


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

%%% Select Frame used as reference and ROI %%%
%%% i_StartForcePeak: [1 3]
%%% i_MaxForcePeak:   [1 6]
%%% i_StartFluoPeak:  [2 3]
%%% i_MaxFluoPeak:    [2 6]
CentralFrame = [1 3];
%%% ROI to use for peaks
% ROISelected  = 4; %[small motor]
ROISelected  = 7; %[whole]
%%% Number of Frame to Take around the CentralFrame
FsFluo = dataGCamp.Info.Fs_FluoImages;
BeforeCentralFrame = 0.32*FsFluo; %(0.32sec -> 8 frames)
AfterCentralFrame  = 2*FsFluo; %(2sec)
%%%
downsamplingfactor = 3; %if FsFluo == 25
% downsamplingfactor = 12; %if FsFluo == 100
%%%

%%% Info images store %%%%%%%%%%%%%%%%%%
UsbPort = 'I';
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


%%% data folder %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
cd('_load_animal_fun')
[folderTASK_FLUO ~] = fun_loadGCaMPanimal(UsbPort,Animal_Name,TrialDay);
cd(CurrDir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Time and Frequencies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t      = dataGCamp.t;
Fs     = dataGCamp.Info.Fs;
FsFluo = dataGCamp.Info.Fs_FluoImages; %[Hz]
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
    Im = imread([folderTASK_FLUO,'\',nameImage]);
    Im_Original = Im;
    
    %filtering
    Im_Original = medfilt2(Im_Original);
    
    %crop image
    Im = Im_Original(...
        ROI_PosAndSize(2,ROISelected):ROI_PosAndSize(2,ROISelected)+ROI_PosAndSize(4,ROISelected),...
        ROI_PosAndSize(1,ROISelected):ROI_PosAndSize(1,ROISelected)+ROI_PosAndSize(3,ROISelected),...
        :);
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
            Im = imread([folderTASK_FLUO,'\',nameImage]);
            Im_Original = Im;
            
            %filtering
            Im_Original = medfilt2(Im_Original);
            
            %crop image
            Im = Im_Original(...
                ROI_PosAndSize(2,ROISelected):ROI_PosAndSize(2,ROISelected)+ROI_PosAndSize(4,ROISelected),...
                ROI_PosAndSize(1,ROISelected):ROI_PosAndSize(1,ROISelected)+ROI_PosAndSize(3,ROISelected),...
                :);
            Im = double(Im);
            
            Im = medfilt2(Im);
            
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
            
            
            %
%             Im_norm(360:end,268:end) = 0;
%             Im_norm(440:end,177:end) = 0;
%             Im_norm(428:end,380:end) = 0;
            
            
        end
        
        %store
        indexStore = indexStore+1;
        MatrixImageSequence_buf(:,:,indexStore) = Im_norm;
        

%             figure
%             imagesc(Im_norm)
%             colormap hot
%             pause(0.5)
%             close
        
    end
    
    MatrixImageSequence{stf,1} = MatrixImageSequence_buf;
    
end

%%%% save %%%%
ImageSequence = save_MatrixImageSequence_fun(dataGCamp.Info.Name,...
                                             dataGCamp.Info.Date,...
                                             MatrixImageSequence,...
                                             CentralFrame,...
                                             ROISelected,...
                                             BeforeCentralFrame,...
                                             AfterCentralFrame,...
                                             downsamplingfactor);
%%%%%%%%%%%%%%

display('End Process')

clear all
 
