%script per estrazione segnale FLUO da diverse ROI e segnale Forza
%il risultato � la struttura dataMouseGCamp_(nome animale)_(giorno test)
%con i segnali e altre info 
%N.B. modificare fun_loadGCaMPanimal per il path

clear
close all
clc


%% Choice of the animal and trial day
UsbPort = 'I';
Animal_Name = 'GCaMP4_control';
TrialDay = '01';
%region of interest to take
ROI_num  = 7;
%size ROIs
ROI_size = [250 250; 250 250; 250 250; 125 125; 125 125; 125 125; 500 500];
ROI_name = {'Motor','Barrel','Occipital','Motor Small','Barrel Small','Occipital Small','Whole'};
%%

%data folder
[folderTASK_FLUO folderTASK_ForceFileName] = fun_loadGCaMPanimal(UsbPort,Animal_Name,TrialDay);
%load force file
dataRobot = importdata(folderTASK_ForceFileName);


%%%%%% data Robot %%%%%%%%%%%%

synchroTask = find(dataRobot(:,7)==0);
%To synchronize robot file (the previous time point)
synchroTask = synchroTask(1)-1;

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
pos_robot = []
if size(dataRobot,2)>8
    %position
    pos_robot = dataRobot(RealStart_Robot:RealStart_Robot+RealDurTask_Robot-1,9);
    %speed
    speed_robot = dataRobot(RealStart_Robot:RealStart_Robot+RealDurTask_Robot-1,10);  
    %acceleration
    acc_robot = dataRobot(RealStart_Robot:RealStart_Robot+RealDurTask_Robot-1,11);
end
%%%%%%


%%%%%% Take task images  %%%%%%

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

%check length of the files
if NumImages<RealDurTask_Fluo
    %difference
    diflen = RealDurTask_Fluo - NumImages;
    %differenza espressa in numero di punti del file del robot
    diflenR = diflen*ResampleParRobotFluo;
    
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
    
    RealDurTask_Fluo = NumImages; 

end


downsamplingfactor=1;

index=0;
wb = waitbar(0,'Images Loading, Please wait...');

%load image stack
% N.B. size of StoreImage depends on num of ROIs
StoreImageAllROIs = cell(ROI_num,1);
ROI_PosAndSize    = zeros(ROI_num,4);

LenFLUO_Im = RealStart_Fluo+RealDurTask_Fluo-1;
for i=RealStart_Fluo :downsamplingfactor:   LenFLUO_Im 
    
    index = index+1;
    indexImage = i+2;
    waitbar(indexImage/(LenFLUO_Im),wb)
    
    nameImage = d(indexImage,1).name;
    
    Im = imread([folderTASK_FLUO,'\',nameImage]);
    Im_Original = Im;
    
    for nR = 1:ROI_num
        
        if index==1
           
            
            %%%%%%%%%%%% region of interest %%%%%%%%%%%%
            fROI = figure('Name','Indicate the ROI','NumberTitle','off');
            hIm = imshow(Im_Original);

            %starting point
              %whole frame
            if strcmp(ROI_name{nR},'Whole')
                hold off
                rect_ROI_default_pos =   [3   3];
            else
                %ROIs
                if nR < 4
                    %first three areas
                    rect_ROI_default_pos =   [116   116]; %-> punto di partenza, area CFA posteriore
                else
                    rect_ROI_default_pos =   [ROI_PosAndSize(nR-3,1)+ROI_PosAndSize(nR-3,3)/4  ROI_PosAndSize(nR-3,2)+ROI_PosAndSize(nR-3,4)/4];  %-> punto centrale previous big square
                    h_rect_Big = imrect(gca, [ROI_PosAndSize(nR-3,1) ROI_PosAndSize(nR-3,2)  ROI_PosAndSize(nR-3,3) ROI_PosAndSize(nR-3,4)]); %-> previous big square
                    hold on
                end
            end
            
            h_rect = imrect(gca, [rect_ROI_default_pos(1) rect_ROI_default_pos(2) ROI_size(nR,1) ROI_size(nR,2)]); %-> finestra di 150x150 o pi� grande px
            pause
            rect_ROI = round(getPosition(h_rect));
            
            
            %         [I_ROI rect_ROI] = imcrop(hIm, rect_ROI);
            %         rect_ROI = round(rect_ROI);
            
            Im = Im_Original(...
                rect_ROI(2):rect_ROI(2)+rect_ROI(4),...
                rect_ROI(1):rect_ROI(1)+rect_ROI(3),...
                :);
            close(fROI)
            %info size ROI
            ROI_PosAndSize(nR,:) = rect_ROI;
            %%%%%%
            
            [nr nc] = size(Im);
            pause(2)
            
%             StoreImageAllROIs{nR,1} = zeros(nr,nc,round(RealDurTask_Fluo/downsamplingfactor));
            
%             StoreImage     = zeros(nr,nc,round(RealDurTask_Fluo/downsamplingfactor));
%             StoreImageNorm = zeros(nr,nc,round(RealDurTask_Fluo/downsamplingfactor));
            
            
        else
%             Im = Im(...
%                 rect_ROI(2):rect_ROI(2)+rect_ROI(4),...
%                 rect_ROI(1):rect_ROI(1)+rect_ROI(3),...
%                 :);
            Im = Im_Original(...
                ROI_PosAndSize(nR,2):ROI_PosAndSize(nR,2)+ROI_PosAndSize(nR,4),...
                ROI_PosAndSize(nR,1):ROI_PosAndSize(nR,1)+ROI_PosAndSize(nR,3),...
                :);
        end
        
        
        %filter (salt and pepper removal)
        Im = medfilt2(Im);
        
%         StoreImage(:,:,index) =  Im;

        StoreImageAllROIs{nR,1}(:,:,index) = Im;
    end
    
            
end
clear d dataRobot Im Im_Original
close(wb)

%%%% choose the baseline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% plot time course  %%%
fbas = figure('Name','Select where animal is not moving');
plot(t_robot,status_robot/10,'r')
hold on
plot(t_robot,fx_robot,'b')

StoreMean     = zeros(ROI_num,LenFLUO_Im);
imFluoResMean = zeros(ROI_num,LenFLUO_Im*ResampleParRobotFluo);

for nR = 1:ROI_num
    %compute average on the areas 
    StoreMean(nR,:)     = squeeze(mean(mean(StoreImageAllROIs{nR,1}(:,:,:),2),1));
    imFluoResMean(nR,:) = resample(StoreMean(nR,:),ResampleParRobotFluo,1);
end
plot(t_robot,imFluoResMean/100);
legend([{'Status','Force'},ROI_name])

[x_1 x_2] = ginput;
%estremi intervallo di interesse
IntervalToFluoMean = round([x_1(1) x_1(2)]*fs_fluo);
close(fbas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%% compute mean baseline Fluo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for nR = 1:ROI_num
    
    StoreImageAllROIs{nR,1} = StoreImageAllROIs{nR,1};
     
    StoreImageMean = round(mean(double(StoreImageAllROIs{nR,1}(:,:,IntervalToFluoMean(1): IntervalToFluoMean(2))),3));
    StoreImageMean(StoreImageMean==0) = 1;
   
    %normalization
    for j=1:size(StoreImageAllROIs{nR,1},3)
        %(f-f0)/f0 
        StoreImageAllROIs{nR,1}(:,:,j) = ((double(StoreImageAllROIs{nR,1}(:,:,j))-StoreImageMean)./StoreImageMean)*100;
        display(j/size(StoreImageAllROIs{nR,1},3)*100)
    end
    
    %average
    m1=mean(StoreImageAllROIs{nR,1},1);
    m2=mean(m1,2);
    %average over the time
    MeanImFluoOverTime = squeeze(m2);
    clear m1 m2    
    
    %resample TimeMeanImage
    imFluoRes = resample(MeanImFluoOverTime,ResampleParRobotFluo,1);
    
    %%%to save data
    % Mean Fluo in ROI i
    dataGCamp.fluoROI(:,nR) = imFluoRes;   
    
    pause(5)
end


%%%%%%%%%%%%%%%
%% PLOT DATA %%
%%%%%%%%%%%%%%%
figure
plot(t_robot,status_robot/10,'r')
hold on
plot(t_robot,fx_robot,'b')
plot(t_robot,dataGCamp.fluoROI/100)
legend([{'Status','Force'},ROI_name])
%%%%


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
    WholeTrials_Index = [1; IndexStatusMENO1([1; find(diff(IndexStatusMENO1)>1)])];
    %l'ultimo valore corrisponde all'ultimo -1 (non re inizia nessuna prova)
    WholeTrials_Index = WholeTrials_Index(1:end-1);
% number of completed trials
    NumTrials = length(WholeTrials_Index) - 1* (status_robot(end)~=-1);
dataGCamp.InfoTrial.NumTrials = NumTrials;


    TrialsVector = zeros(length(status_robot),1);
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
dataGCamp.Info.Name = Animal_Name;
% day
dataGCamp.Info.Date = TrialDay;
% treatment
if strfind(Animal_Name,'stroke')
    dataGCamp.Info.Treat = 1;
elseif strfind(Animal_Name,'control')
    dataGCamp.Info.Treat = 0;
end
% frequency of sampling
dataGCamp.Info.Fs = fs_robot;

% interval Fluo and Force Mean
dataGCamp.Info.IntervalToFluoForceMean = IntervalToFluoMean;

% ROIs
dataGCamp.Info.ROI = ROI_name;
dataGCamp.Info.ROI_PosAndSize = ROI_PosAndSize';
%%%%%%

%%%%%%
filename_GCamp = ['dataMouseGCamp_',Animal_Name,'_',TrialDay ];
save(filename_GCamp,'dataGCamp')
%%%%%%





