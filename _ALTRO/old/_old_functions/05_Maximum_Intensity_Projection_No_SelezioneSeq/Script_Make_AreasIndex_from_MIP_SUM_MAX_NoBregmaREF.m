%%%% MIP %%%
%%% find common area from MIP %%%
clear all
clc

%%% select animal
Animal_name = 'GCaMPChR2_1_control';
%%%

%%%
CurrDir = cd;
Dir_Store_MIP_MAT = '_Store_MIP_all_Days_and_MAT\NoBregmaREF';
MIP_tech_List     = {'max','sumA','sumP'};
LenMIP            = length(MIP_tech_List);

RegionBound_Index = cell(LenMIP, 2);

%%%MAX%%% %%% %%% %%% %%% %%% %%%
if exist([Dir_Store_MIP_MAT,'\',Animal_name] )
    
    for i=1:LenMIP 
        
        %load Im_med
        if strcmp(MIP_tech_List{i}(1:3),'sum')
            MIP_name = MIP_tech_List{i}(1:end-1);
        else
            MIP_name = MIP_tech_List{i};
        end
        
        load([CurrDir,'\',Dir_Store_MIP_MAT,'\',Animal_name,'\Im_AlongDays_',Animal_name,'_MIP_',MIP_name,'.mat'])
        
        figPl = figure;
        imagesc(Im_med)
        colorbar
        %plot rect to select area
        h_rect = imrect(gca,[10 10 50 50]);
        pause
        rect_ROI = round(getPosition(h_rect));
        
        %take thi spart of image
        I_cut = Im_med(rect_ROI(2):rect_ROI(2)+rect_ROI(4),rect_ROI(1):rect_ROI(1)+rect_ROI(3));
        %black image
        Im_buf = Im_med*0;
        Im_buf(rect_ROI(2):rect_ROI(2)+rect_ROI(4),rect_ROI(1):rect_ROI(1)+rect_ROI(3)) = I_cut;
        
        figPl2 = figure;
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
    
    save([Animal_name,'_MIP_RegionBound_Index'],'RegionBound_Index')
    
else
    error([Dir_Store_MIP_MAT,'\',Animal_name,' does not exist'] )
end


