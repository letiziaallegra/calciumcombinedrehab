%%%% MIP %%%
%%% find common area from MIP %%%
clear all
clc


%%%%%%%% Animal_Name %%%%%%%%%%%%%%%%%%%
% ListAnimalTogether = {  'GCaMPChR2_7_control','GCaMPChR2_17_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
%                         'GCaMPChR2_1_control',...
%                         'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke','GCaMPChR2_25_stroke',...
%                         'GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT',...
%                         'GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT','GCaMPChR2_16_stroke_BoNT'};
ListAnimalTogether = {  'GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT',...
                        'GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT','GCaMPChR2_16_stroke_BoNT'};


%%%%%%%%%
CurrDir = cd;
cd .. 
PrevDir = cd;
Dir_Save_IM       = [CurrDir,'\Store_Figures'];
Dir_Save_data     = [CurrDir,'\Store_MAT'];
MIP_tech_List     = {'max','sumA','sumP'};
LenMIP            = length(MIP_tech_List);
%%%%%%%%%


for LATo = 1:length(ListAnimalTogether)

%%% select animal               
Animal_name  = [ListAnimalTogether{LATo}];
%%%


if strcmp(Animal_name,'GCaMPChR2_1_control')
    Dir_Store_MIP_MAT = '_Store_MIP_all_Days_and_MAT\NoBregmaREF';
else
    Dir_Store_MIP_MAT = '_Store_MIP_all_Days_and_MAT';
end

RegionBound_Index = cell(LenMIP, 2);


%%%MAX%%% %%% %%% %%% %%% %%% %%%
% if exist([PrevDir,'\',Dir_Store_MIP_MAT,'\',Animal_name] )
    
    for i=1:LenMIP 
        
        %load Im_med
        if strcmp(MIP_tech_List{i}(1:3),'sum')
            MIP_name = MIP_tech_List{i}(1:end-1);
        else
            MIP_name = MIP_tech_List{i};
        end
        
        if strfind(Animal_name,'stroke_BoNT')
            WeekStr = 'Week_4_';
        else
            WeekStr = [];
        end
        
        load([PrevDir,'\',Dir_Store_MIP_MAT,'\',Animal_name,'\Im_AlongDays_',WeekStr,Animal_name,'_MIP_',MIP_name,'.mat'])
        
        
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
        subplot(321)
        imagesc(Im_buf)
        colorbar
                        
        for MX =1:max(max(Im_buf))
            
            Im_buf_MX = Im_buf>=MX;
        
            subplot(3,2,MX+1)
            imagesc(Im_buf_MX)
            colorbar 
            
            title(num2str(MX))
        
            [bound,label] = bwboundaries(Im_buf_MX);
            temp1 = regionprops(label,'Centroid');
        
            %%% Max Area
            [n in] = max(cellfun('length',bound));
    
            RegionBound_Index{i,MX} = bound{in,1};
                       
        end
        namefig = ['Fig_Thresh_',MIP_tech_List{i},'_MIP_',Animal_name];
        saveas(figPl2,[Dir_Save_IM ,'\',namefig],'fig');
        saveas(figPl2,[Dir_Save_IM ,'\',namefig],'jpeg');
        RegionBound_Index{i,6} = [MIP_tech_List{i},'_MIP'];
       
        close(figPl)
        close(figPl2)
    
    end
    
    save([Dir_Save_data,'\',Animal_name,'_',WeekStr,'MIP_Region_Index'],'RegionBound_Index')
    
% else
%     error([Dir_Store_MIP_MAT,'\',Animal_name,' does not exist'] )
% end

end %end LaTO
