%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Script per trovare le common areas di le MIP/SIP         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

%%% select animal
MainDir = '/Users/alessandro/Desktop/180424_RehabOptogen/MATLAB/';
MainDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB/';
MainDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_GCaMP27-29/';
MainDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_OR/';
MainDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_CTRL_190122/';
MainDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_CTRL_190117/';
% Animal_name = 'OR9_wk01_sham';
[~, Animal_name] = fileparts(uigetdir(MainDir));

%%%



%%%
CurrDir = '';
Dir_Store_MIP_MAT = '/Users/alessandro/Desktop/ELABORAZIONE DATA/05_Maximum_Intensity_Projection_SeqSelez/_Store_MIP_SIP_all_Days_and_MAT/';
Dir_Store_MIP_MAT_RB = '/Users/alessandro/Desktop/ELABORAZIONE DATA/05_Maximum_Intensity_Projection_SeqSelez/_Store_MIP_SIP_RegionBound/';
B00_MakeMIPSIPFolders([MainDir,Animal_name],Dir_Store_MIP_MAT)



% Analyses list - leave it
MIP_tech_List     = {'MIP','SIP_Rostral','SIP_Caudal'};
LenMIP            = length(MIP_tech_List);

RegionBound_Index = cell(LenMIP, 2);

%%%MAX%%% %%% %%% %%% %%% %%% %%%
if exist([Dir_Store_MIP_MAT,filesep,Animal_name], 'dir')
   
    if ~isempty(strfind(Animal_name,'BoNT')) || ~isempty(strfind(Animal_name,'Rehab')) || ~isempty(strfind(Animal_name,'roboti'))
        % list of weeks !!!! IT IS IMPORTANT THAT THE WEEKDAY NUMBER IS
        % MATCHES THE REAL NUMBER, I.E. DAY #6 => FIRST DAY OF SECOND WEEK
        Str_Week = {'Week_1_'; 'Week_2_'; 'Week_3_'};
        % FOR ONE WEEK IT DOES NOT MATTER
        % Str_Week = {'Week_4_'};
        LenWeek = length(Str_Week);
    else
        Str_Week = {''};
        LenWeek = 1;
    end
    
    for iw=1:LenWeek %for week
        
        for i=1:LenMIP           
            
            %load Im_med
            if strcmp(MIP_tech_List{i}(1:3),'SIP')
                MIP_name = MIP_tech_List{i}(1:3);
            else
                MIP_name = MIP_tech_List{i};
            end
            
            %load Im_med
            load([CurrDir,filesep,Dir_Store_MIP_MAT,filesep,Animal_name,filesep,'Im_AlongDays_', Str_Week{iw}, Animal_name,'_Long_SelSeq_', MIP_name,'.mat'])
            
            figPl = figure('Name',['Cut the rectangular area out (', MIP_tech_List{i},')' ]);
            imagesc(Im_med)
            colorbar
            %plot rect to select area
            switch MIP_tech_List{i}
                case 'MIP'
                    h_rect = imrect(gca,[1 1 size(Im_med,1)-1 size(Im_med,2)-1]);
                case 'SIP_Rostral'
                    h_rect = imrect(gca,[1 1 size(Im_med,1)-1 size(Im_med,2)/2]);
                case 'SIP_Caudal'
                    h_rect = imrect(gca,[1 size(Im_med,2)/2 size(Im_med,1)-1 size(Im_med,2)-size(Im_med,2)/2]);
                case 'SIP'
                    h_rect = imrect(gca,[1 1 size(Im_med,1)-1 size(Im_med,2)/2]);
            end
            
            pause
            rect_ROI = round(getPosition(h_rect));
            
            %take this part of image
            I_cut = Im_med(rect_ROI(2):rect_ROI(2)+rect_ROI(4),rect_ROI(1):rect_ROI(1)+rect_ROI(3));
            %black image
            Im_buf = Im_med*0;
            Im_buf(rect_ROI(2):rect_ROI(2)+rect_ROI(4),rect_ROI(1):rect_ROI(1)+rect_ROI(3)) = I_cut;
            
            figPl2 = figure('Name',['Select ', MIP_tech_List{i}]);
            subplot(121)
            imagesc(Im_buf)
            colorbar
            
            [x y] = ginput;
            MX = Im_buf(round(y),round(x));
            Im_buf_MX = Im_buf>=MX;
            
            subplot(122)
            imagesc(Im_buf_MX)
            colorbar
            pause(1)
            
            [bound,label] = bwboundaries(Im_buf_MX);
            temp1 = regionprops(label,'Centroid');
            
            %%% Max Area
            [n in] = max(cellfun('length',bound));
            
            RegionBound_Index{i,1} = bound{in,1};
            RegionBound_Index{i,2} = [MIP_tech_List{i},'_MIP'];
            
            close(figPl)
            close(figPl2)
            
        end
        cdir = pwd;
        %save RegionBound_Index
        st_dir = [CurrDir,filesep,Dir_Store_MIP_MAT_RB,filesep,Animal_name,filesep];
        if ~exist(st_dir,'dir')
            mkdir(st_dir)
        end
        cd(st_dir)
        save([Animal_name,'_MIP_SIP_',Str_Week{iw},'RegionBound_Index'],'RegionBound_Index')
        
        cd(cdir)
    end%for week
    
    
    
else
    error([Dir_Store_MIP_MAT,filesep,Animal_name,' does not exist'] )
end

%OUTPUT MOVE TO _Store_MIP_SIP_RegionBound/AnimalName
