%
% script to syncrhonize fluo and force signals and
% define the ROIs signals based on the Centroid and the Radius assessed
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Radius assessed
Radius        = round(167.51/2); %[pixels] <- FISSO
%Centroids
Centroid_List = round([340 418; 399 380]); %[pixels] <- MODIFICABILE (controlla Centroid_of_ROI_Animals)
ROI_num_Centroid = size(Centroid_List,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Info images store %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UsbPort = 'H';
%
Animal_Name = dataGCamp.Info.Name;
TrialDay    = dataGCamp.Info.Date;
%interval of frames to compute f0
IntervalToFluoMean = dataGCamp.Info.IntervalToFluoForceMean_PointFluoFreq; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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




% %storage of each f0 for each ROI
% MatrixImage_f0_ROI = zeros(sz_Im(1),sz_Im(2),ROI_num_Centroid);

index=0;
downsamplingfactor=1;

wb = waitbar(0,'Images Loading, Please wait...');

%%% load data Images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LenFLUO_Im = RealStart_Fluo+RealDurTask_Fluo-1;
for i=RealStart_Fluo :downsamplingfactor:   LenFLUO_Im 
    
    index = index+1;
    indexImage = i+2;
    waitbar(indexImage/(LenFLUO_Im),wb)
    
    nameImage = d(indexImage,1).name;
    
    Im = imread([folderTASK_FLUO,'\',nameImage]);
    Im_Original = Im;
    sz_Im = size(Im_Original);
    
    %filtering
    Im_Original = medfilt2(Im_Original);  
    
    for nR = 1:ROI_num_Centroid
        
        if index==1           
                       
            %%%%%%%%%%%% region of interest based on Centroid and Radius %%%%%%%%%%%%
            fROI = figure('Name','Show the ROIs','NumberTitle','off');
            hIm = imshow(Im_Original);
            
            
            % Make a logical image with the selected circular region set to 1, the rest % to zero
            [xgrid, ygrid] = meshgrid(1:sz_Im (2), 1:sz_Im (1));
            
            Cx = xgrid - Centroid_List(nR,1);    % offset the origin
            Cy = ygrid - Centroid_List(nR,2);
            circlemask = Cx.^2 + Cy.^2 <= Radius.^2;
            
            
            %%% Trova f0 di questa area circolare %%%
            MatrixImage_f0 = zeros(sz_Im(1),sz_Im(2),IntervalToFluoMean(2)-IntervalToFluoMean(1)+1);
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
                circle_image_f0 = double(Im_Original_f0) .* circlemask;
                circle_image_f0(circle_image_f0==0) = NaN;
                
                %
                indexStore_f0 = indexStore_f0+1;
                MatrixImage_f0(:,:,indexStore_f0) = circle_image_f0;
                
            end
            display('F0 extracted')
            
           
            MeanMatrix_f0 = nanmean( double(MatrixImage_f0),3); %in gray tones
            MeanMatrix_f0(MeanMatrix_f0==0) = 1;
            
            figure
            imshow(MeanMatrix_f0, []);
            
            %storage delle f0 delle ROI
            MatrixImage_f0_ROI(:,:,nR) = MeanMatrix_f0;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
            
        end
            
            % Use the mask to select part of the image
            circle_image = double(Im_Original) .* circlemask;                       
            %Metto tutti gli 0 pari a NaN
            circle_image(circle_image==0) = NaN;
            
            %fai delta_f/f0
            Im = (circle_image-squeeze(MatrixImage_f0_ROI(:,:,nR)))./squeeze(MatrixImage_f0_ROI(:,:,nR))*100;            
           
            %compute average of the image (just the circle bc the outside the circle is all NaN
            ROI_Centroid_MeanSignal(i,nR)     = squeeze(nanmean(nanmean(Im,2),1));
            
    end
end

imFluoResMean = [];
for nR = 1:ROI_num_Centroid
    
    imFluoResMean_Buf = resample(ROI_Centroid_MeanSignal(:,nR),ResampleParRobotFluo,1);
    imFluoResMean_Buf = imFluoResMean_Buf - median(imFluoResMean_Buf);
    
    imFluoResMean     = [imFluoResMean imFluoResMean_Buf];
end
  

%%%%%%%%%%%%%%%
%% PLOT DATA %%
%%%%%%%%%%%%%%%
figure
plot(t,dataGCamp.status/10,'r')
hold on
plot(t,dataGCamp.fx,'b')
plot(t,imFluoResMean/100)
legend([{'Status','Force'}])
%%%%


%%%%%%%%%%%%%%%
%% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%

% image f0 ROI
dataGCamp.CentroidAnalysis.ImageROI_f0   = MatrixImage_f0;
% centroid List
dataGCamp.CentroidAnalysis.Centroid_List = Centroid_List;
% Radius
dataGCamp.CentroidAnalysis.Radius        = Radius;
% signals
dataGCamp.CentroidAnalysis.ROI_Signal    = imFluoResMean;


%%%%%%
filename_GCamp = ['dataMouseGCamp_',dataGCamp.Info.Name,'_',dataGCamp.Info.Date,'_Par','_Centroid' ];
save(filename_GCamp,'dataGCamp')
%%%%%%





