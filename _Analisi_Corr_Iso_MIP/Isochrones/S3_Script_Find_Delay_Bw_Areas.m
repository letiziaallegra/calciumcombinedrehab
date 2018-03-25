%Delays among areas (isochrones) -> plot

clear
% close all
clc


HD = 'F';
User = getenv('USERNAME');

%MEAN or MEDIAN
% MeaMed = 1; %MEAN
MeaMed = 2; %MEDIAN


%%%%%%%%%%%%%%%%
if MeaMed == 1
    ISO_FOLDER  = [HD,':\LENS\Isocrone_SelSeq\Figures_Weekly_Mean of Trials\'];
elseif MeaMed == 2
    ISO_FOLDER  = [HD,':\LENS\Isocrone_SelSeq\Figures_Weekly_Median of Trials\'];
end

TREAT       = {'control','stroke','rehab'};
TREAT_NAME  = {'control','stroke','rehab','rehab4w'};
%%%%%%%%%%%%%%%%

%%%%%%%% load the list of indeces corresponding to the functional areas%
load('RegionArea_FunctionalMask_Index'); %-> RegionArea_FunctionalMask_Index
AREA_NAME   = {'CFA','M1Sh','M1HL','M1Tk','S1','S1Bc','RSD','V1'};
FuncAres_ToTake      = {2, 3, 4, 5, 6, 7, 8, 10};
% AREA_NAME   = {'M1 CFA','otherM1','S1','S1Bc','RSD','V1'};
% FuncAres_ToTake      = { 2, [3, 4, 5], 6, 7, 8, 10};
Num_FunAreas  = size(RegionArea_FunctionalMask_Index,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% list animal for plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_N_control   = 'GCaMPChR2_7_control_Week_1';
P_N_stroke    = 'GCaMPChR2_25_stroke_Week_1';
P_N_rehab_1W  = 'GCaMPChR2_14_stroke_BoNT_Week_1';
P_N_rehab_4W  = 'GCaMPChR2_14_stroke_BoNT_Week_4';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DD_Matrix = cell(4,1);

for i_t = 1:length(TREAT) %for treat
    
    F = [ISO_FOLDER,'\',TREAT{i_t}];
    file_F = dir(F);
    
    DD_Matrix_Treat    = [];
    DD_Matrix_Treat_w4 = [];
    
    for i_f=3:length(file_F) %for file
         
        fileAnimal = file_F(i_f,1).name;
        W1 = ~isempty(strfind(fileAnimal,'Week_1'));
        W4 = ~isempty(strfind(fileAnimal,'Week_4'));
        
        
        if ((((i_t == 1) | (i_t == 2)) && ~isempty(strfind(fileAnimal,'GCaMP')) && ~isempty(strfind(fileAnimal,'.fig')))  | ...  %if GCAMP .fig
           ( (i_t == 3) && (W1 | W4)  && ~isempty(strfind(fileAnimal,'GCaMP')) && ~isempty(strfind(fileAnimal,'.fig'))) ) &... 
           ( isempty(strfind(fileAnimal,'GCaMPChR2_16_stroke_BoNT')) )%if GCAMP 
%        
%        if ((((i_t == 1) | (i_t == 2)) && ~isempty(strfind(fileAnimal,'GCaMP')) && ~isempty(strfind(fileAnimal,'.fig')))  | ...  %if GCAMP .fig
%            ( (i_t == 3) && (W1 | W4)  && ~isempty(strfind(fileAnimal,'GCaMP')) && ~isempty(strfind(fileAnimal,'.fig'))) ) &... 
%            ( isempty(strfind(fileAnimal,'GCaMP16_stroke_BoNT')) &&  isempty(strfind(fileAnimal,'GCaMP17_stroke_BoNT')) &&  isempty(strfind(fileAnimal,'GCaMP18_stroke_BoNT'))) &&  isempty(strfind(fileAnimal,'GCaMP16_ChR2_16_stroke_BoNT'))%if GCAMP 
%        
%         if (((i_t == 1) | (i_t == 2)) && ~isempty(strfind(fileAnimal,'GCaMP')) && ~isempty(strfind(fileAnimal,'.fig')))  | ...  %if GCAMP .fig
%            ( (i_t == 3) && (W1 | W4)  && ~isempty(strfind(fileAnimal,'GCaMP')) && ~isempty(strfind(fileAnimal,'.fig')))  
% %             
            %%%load file
            %load([F,'\',file_F]) %-> Fig Matlab
            open([F,'\',fileAnimal]);
            h = gca;
            M_Iso = getimage(h);
            close
            
            %to return to the 32x32 areas (set from Script_for_Isochrone_Weekly)
            M_Iso_32 = imresize(M_Iso,1/40);
            %to return to the original resolution 512x512
            M = imresize(M_Iso_32,512/32);
            
            [s1 s2 s3] = size(M);
            
            Matrix_Delay_An = zeros(length(FuncAres_ToTake),1);
            
            %%% Find Areas Delay for the animal %%%
            for i_fAtT = 1:length(FuncAres_ToTake)
                Area_index = FuncAres_ToTake{i_fAtT};
                if length(Area_index) == 1
                    M_Areas = M(sub2ind([s1 s2],RegionArea_FunctionalMask_Index{Area_index ,1}(:,2),RegionArea_FunctionalMask_Index{Area_index ,1}(:,2)));
                    
                    if MeaMed == 1                        
                        Matrix_Delay_An(i_fAtT,1) = nanmean(M_Areas);
                    elseif MeaMed == 2
                        Matrix_Delay_An(i_fAtT,1) = nanmedian(M_Areas);
                    end
                    
%                     %%%plot test
%                     M_Plot = M;
%                     M_Plot(sub2ind([s1 s2],RegionArea_FunctionalMask_Index{Area_index ,1}(:,2),RegionArea_FunctionalMask_Index{Area_index ,1}(:,1))) = 10000;
%                     imagesc(M_Plot)
%                     ciao = 1;
%                     %%%  
                    
                else
                    %area composed by more than one simple area
                    M_Areas = [];
                    for  i_AInd = 1:length(Area_index)                        
                        M_Areas = [M_Areas; (M(sub2ind([s1 s2],RegionArea_FunctionalMask_Index{Area_index(i_AInd) ,1}(:,1),RegionArea_FunctionalMask_Index{Area_index(i_AInd)  ,1}(:,2))))];
                    end
                    
                    if MeaMed == 1
                        Matrix_Delay_An(i_fAtT,1) = nanmean(M_Areas);
                    elseif MeaMed == 2
                        Matrix_Delay_An(i_fAtT,1) = nanmedian(M_Areas);
                    end
                end
            end
            %%%
            
%             %%% difference bw M1 and other areas
%             M1_time =  Matrix_Delay_An(1);
%             Matrix_Delay_An  = (Matrix_Delay_An - M1_time);
%             %             Matrix_Delay_An  = abs(Matrix_Delay_An-M1_time);
%             %%%
%             
%             %%% store delay diff
%             if W1
%                 %%%
%                 DD_Matrix_Treat = [DD_Matrix_Treat, [Matrix_Delay_An(2:end)]];
%             elseif W4
%                 %%%
%                 DD_Matrix_Treat_w4 = [DD_Matrix_Treat_w4, [Matrix_Delay_An(2:end)]];
%             end
%             %%%
                        
             %%% tutti  %%%
            Matrix_Delay_An  = (Matrix_Delay_An);
            %%% store delay diff
            if W1
                %%%
                DD_Matrix_Treat = [DD_Matrix_Treat, [Matrix_Delay_An(1:end)]];
            elseif W4
                %%%
                DD_Matrix_Treat_w4 = [DD_Matrix_Treat_w4, [Matrix_Delay_An(1:end)]];
            end
            %%%
            
            
            %save for PLOT
            if strfind(fileAnimal,P_N_control)                
                M_Iso_Store(:,:,1)  = M;                
            elseif strfind(fileAnimal,P_N_stroke)                
                M_Iso_Store(:,:,2)  = M;                
            elseif strfind(fileAnimal,P_N_rehab_1W)
                M_Iso_Store(:,:,3)  = M;
            elseif strfind(fileAnimal,P_N_rehab_4W) 
                M_Iso_Store(:,:,4)  = M;
            end
            
            
            
        end
        
    end %end for file
    
    
    DD_Matrix{i_t,1} = DD_Matrix_Treat;    
    
    if (i_t == 3) && (W4)
        DD_Matrix{4,1} = DD_Matrix_Treat_w4;
    end
        
    
end


%%%plot%%%
COLOR_T =  colormap(hsv(4));
figure
subplot(1,5,1)
% subplot(2,4,[5:8])
for i_T=1:4
%     
    Med = nanmean(DD_Matrix{i_T,1},2);
%     Med = nanmedian(DD_Matrix{i_T,1},2);
    SDE = nanstd(DD_Matrix{i_T,1},[],2)/sqrt(size(DD_Matrix{i_T,1},2));
    errorbar(Med,SDE,'color',COLOR_T(i_T,:)) 

% subplot(2,2,i_T)
% plot(DD_Matrix{i_T,1})    

    hold on
    
    title(['Delay'])
%     xlim([0.5 11.5])
    set(gca,'xTick',     1:length(AREA_NAME(2:end)))
    set(gca,'xTickLabel',AREA_NAME(2:end))
%     ylim([0.1 1.1])
%     set(gca,'yTick',     0.1:0.1:1)
    ylabel('Delay_{Area} - Delay_{M1}')
    xlabel('Area')
    
end
rotateXLabels( gca, 60 )
legend(TREAT_NAME(1:4),'Location','Best')


% PosPlot              = [2 3 4 5]
% % PosPlot              = [1:4]
% for i_T=1:4
%     
%     M_Iso   = M_Iso_Store(:,:,i_T);
%     for i_funA=1:length(FuncAres_ToTake)
%         
%         M_Iso_Buf = zeros(size(M_Iso_Store(:,:,i_T),1),size(M_Iso_Store(:,:,i_T),2));
%         
%         Area_index = FuncAres_ToTake{i_funA};   
%         for  i_AInd = 1:length(Area_index)
%             M_Iso_Buf(sub2ind([s1 s2],RegionArea_FunctionalMask_Index{Area_index(i_AInd),1}(:,2),RegionArea_FunctionalMask_Index{Area_index(i_AInd),1}(:,1))) = 1;
%         end
%         %%% Max Area
%         M_Iso_Buf=  imfill(M_Iso_Buf,'holes');
%         se = strel('disk',4);        
%         M_Iso_Buf = imdilate(M_Iso_Buf,se);
%         [bound,label] = bwboundaries(M_Iso_Buf);
%         [n in] = max(cellfun('length',bound));
%         M_Iso(sub2ind([s1 s2],bound{in,1}(:,1),bound{in,1}(:,2))) = 0.5;
%     end
%     M_Iso_PLOT = M_Iso;
%     
% %     %to return to the 32x32 areas (set from Script_for_Isochrone_Weekly)
% %     M_Iso_1080 = imresize(M_Iso,40);
% %     %to return to the original resolution 512x512
% %     M_Iso_PLOT = imresize(M_Iso_1080,32/512);
%     
%     M_Iso_PLOT(isnan(M_Iso)) = -0.5;
%     se = strel('disk',3);        
%     M_Iso_PLOT = imdilate(M_Iso_PLOT,se);
%     M_Iso_PLOT( M_Iso_PLOT==100) = NaN;
%     
% %     subplot(2,2,PosPlot(i_T))
%     subplot(1,5,PosPlot(i_T))
%     imagesc(M_Iso_PLOT)
%     
%     changeLabel_MIP(size(M_Iso_PLOT,1),size(M_Iso_PLOT,2),1);
%     caxis([-0.1 0.3])
%     
%     title(TREAT_NAME(i_T));
% 
% end
% colorbar
