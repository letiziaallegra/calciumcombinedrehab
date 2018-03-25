%
% script to extract fluo frames corresponding to buzz sound
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
soundStatus = 7;
% soundStatus = 4;
%%% ROI to use for peaks
ROISelected  = 7; %[whole]
%%% Number of Frame to Take around the CentralFrame
CentralFrame       = [];
BeforeCentralFrame = [];
AfterCentralFrame  = [];
%%%
downsamplingfactor = 1; %if FsFluo == 25
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


%%% Time and Status and Frequencies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t      = dataGCamp.t;
status = dataGCamp.status;
Fs     = dataGCamp.Info.Fs;
FsFluo = dataGCamp.Info.Fs_FluoImages; %[Hz]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Find Index of the Buzz Sound %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
status_buzz = find(status==soundStatus);
status_buzz_diff_end   = [status_buzz(diff(status_buzz)>1); status_buzz(end)];
status_buzz_diff_start = [status_buzz(1); status_buzz(circshift(diff(status_buzz)>1,1))];
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
TimeTrigger = [status_buzz_diff_start status_buzz_diff_end]/Fs; %[100 Hz]
StoreTrigger_Fluo = [];
for tt=1:length(TimeTrigger)
    
    Start_Int = round(TimeTrigger(tt,1)*FsFluo);
    End_Int   = round(TimeTrigger(tt,2)*FsFluo);
    
    StoreTrigger_Fluo = [StoreTrigger_Fluo; [Start_Int End_Int]];
    
end
minDimAud = min(abs(StoreTrigger_Fluo(:,2)-StoreTrigger_Fluo(:,1)));
StoreTrigger_Fluo(:,2) = StoreTrigger_Fluo(:,1) +  minDimAud;

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
            
                    %store
        indexStore = indexStore+1;
        MatrixImageSequence_buf(:,:,indexStore) = Im_norm;
        end
        

        

%             figure
%             imagesc(Im_norm)
%             colormap hot
%             pause(0.5)
%             close
        
    end
    
    MatrixImageSequence{stf,1} = MatrixImageSequence_buf;
    
end

%%%% save %%%%
if soundStatus == 2 || soundStatus == 4
ImageSequence = save_MatrixImageSequence_fun_Audio(dataGCamp.Info.Name,...
                                                 dataGCamp.Info.Date,...
                                                 MatrixImageSequence,...
                                                 CentralFrame,...
                                                 ROISelected,...
                                                 BeforeCentralFrame,...
                                                 AfterCentralFrame,...
                                                 downsamplingfactor);
elseif soundStatus == 7
    ImageSequence = save_MatrixImageSequence_fun_Rest(dataGCamp.Info.Name,...
                                                 dataGCamp.Info.Date,...
                                                 MatrixImageSequence,...
                                                 CentralFrame,...
                                                 ROISelected,...
                                                 BeforeCentralFrame,...
                                                 AfterCentralFrame,...
                                                 downsamplingfactor);
end
    
%%%%%%%%%%%%%%

display('End Process')

clear all
 
