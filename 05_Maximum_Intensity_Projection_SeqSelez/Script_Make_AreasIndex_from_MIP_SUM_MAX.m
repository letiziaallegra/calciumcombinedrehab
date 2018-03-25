%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Script per trovare le common areas di le MIP/SIP         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

%%% select animal
Animal_name = 'GCampChR2_TOX5';
%%%

%%%
CurrDir = cd;
Dir_Store_MIP_MAT = '_Store_MIP_SIP_all_Days_and_MAT';
MIP_tech_List     = {'MIP','SIP_Rostral','SIP_Caudal'};
LenMIP            = length(MIP_tech_List);

RegionBound_Index = cell(LenMIP, 2);

%%%MAX%%% %%% %%% %%% %%% %%% %%%
if exist([Dir_Store_MIP_MAT,'\',Animal_name] )
    
    if ~isempty(strfind(Animal_name,'BoNT')) | strfind(Animal_name,'Rehab')
%         Str_Week = {'Week_1_'; 'Week_2_'; 'Week_3_'; 'Week_4_'};
        Str_Week = {'Week_4_'};
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
            load([CurrDir,'\',Dir_Store_MIP_MAT,'\',Animal_name,'\Im_AlongDays_', Str_Week{iw}, Animal_name,'_Long_SelSeq_', MIP_name,'.mat'])
            
            figPl = figure('Name',['Cut the rectangular area out (', MIP_tech_List{i},')' ]);
            imagesc(Im_med)
            colorbar
            %plot rect to select area
            h_rect = imrect(gca,[10 10 50 50]);
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
        
        %save RegionBound_Index
        save([Animal_name,'_MIP_SIP_',Str_Week{iw},'RegionBound_Index'],'RegionBound_Index')
        
    end%for week
    
    
    
else
    error([Dir_Store_MIP_MAT,'\',Animal_name,' does not exist'] )
end


