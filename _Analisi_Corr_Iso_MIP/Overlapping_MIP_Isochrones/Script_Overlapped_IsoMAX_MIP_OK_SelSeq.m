%
% Find percentage of M1 and other areas that constitute the ROI of the MIP
%


clear
close all
clc

CurrDir = cd;


%%%% Choice of the animal and trial day
UsbPortHD = 'F';
UserName  = getenv('username');
AnimalDir = [UsbPortHD,':\LENS\Animals Data'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Animal_Name %%%%%%%%%%%%%%%%%%%
% ListAnimalTogether = { 'GCaMPChR2_7_control'};


% ListAnimalTogether = {  'GCaMPChR2_7_control','GCaMPChR2_17_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
%                         'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke','GCaMPChR2_25_stroke',...   
%                         'GCaMP16_stroke_BoNT', 'GCaMP17_stroke_BoNT', 'GCaMP18_stroke_BoNT'...
%                         'GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT','GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT','GCaMPChR2_16_stroke_BoNT'};
                    
%  ListAnimalTogether = {  'GCaMPChR2_7_control','GCaMPChR2_17_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
%                         'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke','GCaMPChR2_25_stroke',...
%                         ...
%                         'GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT','GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT'};

ListAnimalTogether = {  'GCaMPChR2_7_control','GCaMPChR2_17_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
                        'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke','GCaMPChR2_25_stroke','GCaMPChR2_26_stroke',...   
                        'GCaMP16_stroke_BoNT', 'GCaMP18_stroke_BoNT'...
                        'GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT','GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT'};
                 


%%%%%%%%%%%%%%% Threshold to use %%%%%%%%%
ThreshMIP = 3; %[minimum number of days when MIP/SIP is present]
ThreshIso = []; %[sec]


%%%%%%%%%%%%%%%% ISOCRONE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ISO_FOLDER  = [UsbPortHD,':\LENS\Isocrone_SelSeq_MAX\Figures_Weekly_Mean of Trials'];
TREAT       = {'control','stroke','rehab'};


% %%%%%%%% load the list of indeces corresponding to the functional areas%
% load('RegionArea_FunctionalMask_Index'); %-> RegionArea_FunctionalMask_Index
% % AREA_NAME   = {'M2','CFA','M1Sh','M1HL','M1Tk','S1','S1Bc','RSD','LPTa','V1','AuD'};
% % % FuncAres_ToTake      = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
% AREA_NAME   = {'CFA','other M1','Sensorial','Caudal'};
% FuncAres_ToTake      = {  2, [1, 3, 4, 5], [6, 7], [8, 9, 10]};
% Num_FunAreas  = size(RegionArea_FunctionalMask_Index,1);


% %%%%%%% list animal for plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P_N_control   = 'GCaMPChR2_7_control_Week_1';
% P_N_stroke    = 'GCaMPChR2_25_stroke_Week_1';
% P_N_rehab_1W  = 'GCaMPChR2_14_stroke_BoNT_Week_1';
% P_N_rehab_4W  = 'GCaMPChR2_14_stroke_BoNT_Week_4';


%Store cells
Area_Over_MIP_STORE     = cell(4,1);
Area_Over_SIP_STORE     = cell(4,1);
DistCentroid_MIP_STORE  = cell(4,1);
DistCentroid_SIP_STORE  = cell(4,1);





for LATo = 1:length(ListAnimalTogether)
    
    Animal_Name                    = [ListAnimalTogether{LATo}];
    
    
    %%%%% LOAD data REGION BOUND MIP %%%%%%%%%%%%%%%%%%%%
    ROI_Folder   = ['C:\Users\',UserName,'\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\05_Maximum_Intensity_Projection_SeqSelez\_Store_MIP_SIP_all_Days_and_MAT\',Animal_Name];
        
    Area_Over_MIP       = [];
    Area_Over_SIP       = [];    
    DistCentroid_MIP    = [];
    DistCentroid_SIP    = []; 

    %%%%%%%%%%
    if strfind(Animal_Name,'BoNT')
        
        
        Str_MIP_Week(1) = {'Week_1_'};
        Str_MIP_Week(2) = {'Week_4_'};
        
        Str_Iso_Week(1) = {'Week_1_'};
        Str_Iso_Week(2) = {'Week_4_'};
        
        Num_week  = 2;   
        
    else 
        
        Str_MIP_Week(1) = {[]};
        
        Str_Iso_Week(1) = {'Week_1_'};    
        
        Num_week  = 1; 
                
    end
       
    
    for i_week=1:Num_week %for week
        
        
        %%%%%%%%% MIP and SIP %%%%%%%%%%
        load([ROI_Folder,'\','Im_AlongDays_',Str_MIP_Week{i_week},Animal_Name,'_Long_SelSeq_MIP']);
        Im_MIP    = Im_med;
        Im_MIP_Th = Im_MIP>=ThreshMIP;
        clear Im_med
        %%Centroid
        [bound_Im label_Im] = bwboundaries(Im_MIP_Th);
        Cent_Im             = regionprops(label_Im,'Centroid');
        [n in]              = max(cellfun('length',bound_Im));
        Cent_Im_MIP         = Cent_Im(in,1);
        if isempty(Cent_Im_MIP)
            Cent_Im_MIP = [];
            Cent_Im_MIP.Centroid = [NaN NaN];
        end
        %%%
        
        load([ROI_Folder,'\','Im_AlongDays_',Str_MIP_Week{i_week},Animal_Name,'_Long_SelSeq_SIP']);
        Im_SIP    = Im_med;
        Im_SIP_Th = Im_SIP>=ThreshMIP;
        clear Im_med
        %%Centroid
        [bound_Im label_Im] = bwboundaries(Im_SIP_Th);
        Cent_Im             = regionprops(label_Im,'Centroid');
        [n in]              = max(cellfun('length',bound_Im));
        Cent_Im_SIP         = Cent_Im(in,1);
        if isempty(Cent_Im_SIP)
            Cent_Im_SIP = [];
            Cent_Im_SIP.Centroid = [NaN NaN];
        end
        %%%
        
        
        %%%%%%%%% IsoChrones %%%%%%%%%%
        if strfind(Animal_Name,TREAT{1})
            %control
            open([ISO_FOLDER,'\',TREAT{1},'\',Animal_Name,'_',Str_Iso_Week{i_week},'5_days_average.fig'])
            
%             ThreshIso = 0.1;
            
        elseif strfind(Animal_Name,TREAT{2})
            
            if strfind(Animal_Name,'BoNT')
                %rehab
                open([ISO_FOLDER,'\',TREAT{3},'\',Animal_Name,'_',Str_Iso_Week{i_week},'5_days_average.fig'])
                
%                 ThreshIso = 0.1;
                
            elseif strfind(Animal_Name,TREAT{2})
                %stroke
                open([ISO_FOLDER,'\',TREAT{2},'\',Animal_Name,'_',Str_Iso_Week{i_week},'5_days_average.fig'])
                
%                 ThreshIso = 0.2;
            end
            
        end
        
        h = gca;
        M_Iso = getimage(h);
        close
        %to return to the 32x32 areas (set from Script_for_Isochrone_Weekly)
        M_Iso_32 = imresize(M_Iso,1/40);
        %to return to the original resolution 512x512
        M_Iso_512 = imresize(M_Iso_32,512/32);
        %threshold Iso
        if isempty(ThreshIso)
            meIso  = nanmean(nanmean(M_Iso_512));
            stdIso = 2*nanstd(nanstd(M_Iso_512));
            M_Iso_512_low  = ones(size(M_Iso_512,1),size(M_Iso_512,2));
            M_Iso_512_high = ones(size(M_Iso_512,1),size(M_Iso_512,2));
            M_Iso_512_low(M_Iso_512<=meIso-stdIso)  = 0;
            M_Iso_512_low(isnan(M_Iso_512))         = 0;
            M_Iso_512_high(M_Iso_512>=meIso+stdIso) = 0;
            M_Iso_512_high(isnan(M_Iso_512))        = 0;
            M_Iso_512_T = M_Iso_512_low+M_Iso_512_high;
            M_Iso_Th = M_Iso_512_T==2;
        else
            M_Iso_Th = M_Iso_512<=ThreshIso;
        end
        
        %%Centroid
        [bound_Im label_Im] = bwboundaries(M_Iso_Th);
        Cent_Im             = regionprops(label_Im,'Centroid');
        [n in]              = max(cellfun('length',bound_Im));
        Cent_Im_Iso         = Cent_Im(in,1);
        
        if isempty(Cent_Im_Iso)
            Iso_thresh_ok = 1;
            ThreshIso_Up  = 0.01;
            while Iso_thresh_ok                
                M_Iso_Th = M_Iso_512<=ThreshIso_Up;
                
                %%Centroid
                [bound_Im label_Im] = bwboundaries(M_Iso_Th);
                Cent_Im             = regionprops(label_Im,'Centroid');
                [n in]              = max(cellfun('length',bound_Im));
                Cent_Im_Iso         = Cent_Im(in,1);
                
                if ~isempty(Cent_Im_Iso)
                    Iso_thresh_ok = 0;
                end
                ThreshIso_Up = ThreshIso_Up + 0.01;
            end
            
        end
        %%%     
        
        
        
        %%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% Overlapping %%%
        Over_MIP_Iso       = Im_MIP_Th + M_Iso_Th;
        Tot_Over_MIP       = sum(sum(Im_MIP_Th));
        Tot_Over_Iso       = sum(sum(M_Iso_Th));
        Tot_Over_MIP_Iso   = sum(sum(Over_MIP_Iso))/2;
        Over_MIP_Iso_Union = Over_MIP_Iso>=1;
        Over_MIP_Iso_Inter = Over_MIP_Iso>=2;
%         Area_Over_MIP      = sum(sum(Over_MIP_Iso_Inter));
%         Area_Over_MIP      = sum(sum(Over_MIP_Iso_Inter))/Tot_Over_MIP_Iso;
%         Area_Over_MIP      = sum(sum(Over_MIP_Iso_Inter))/Tot_Over_MIP;
        Area_Over_MIP      = sum(sum(Over_MIP_Iso_Inter))/sum(sum(Over_MIP_Iso_Union));
%         Area_Over_MIP      = sum(sum(Over_MIP_Iso_Inter))/ (Tot_Over_MIP+Tot_Over_Iso);


        if Area_Over_MIP == 0
            Area_Over_MIP = NaN;
        end
        
        Over_SIP_Iso       = Im_SIP_Th + M_Iso_Th;
        Tot_Over_SIP       = sum(sum(Im_SIP_Th));
        Tot_Over_Iso       = sum(sum(M_Iso_Th));
        Tot_Over_SIP_Iso   = sum(sum(Over_SIP_Iso))/2;
        Over_SIP_Iso_Union = Over_SIP_Iso>=1;
        Over_SIP_Iso_Inter = Over_SIP_Iso>=2;
        Area_Over_SIP      = sum(sum(Over_SIP_Iso_Inter));
%         Area_Over_SIP      = sum(sum(Over_SIP_Iso_Inter))/Tot_Over_SIP_Iso;
%         Area_Over_SIP      = sum(sum(Over_SIP_Iso_Inter))/Tot_Over_SIP;
        Area_Over_SIP      = sum(sum(Over_SIP_Iso_Inter))/sum(sum(Over_SIP_Iso_Union));
%         Area_Over_SIP      = sum(sum(Over_SIP_Iso_Inter))/(Tot_Over_SIP+Tot_Over_Iso);

%         figure
%         subplot(221)
%         imagesc(Im_SIP_Th)
%         subplot(222)
%         imagesc(M_Iso_Th)
%         subplot(223)
%         imagesc(Over_SIP_Iso_Union)
%         subplot(224)
%         imagesc(Over_SIP_Iso_Inter)
        
%         if Area_Over_SIP == 0
%             Area_Over_SIP = NaN;
%         end
       
        
%         figure('Name',Animal_Name)
%         subplot(121)
%         imagesc(Over_SIP_Iso_Inter)
%         subplot(122)
%         imagesc(Over_SIP_Iso_Union)
        
        %%%%%%%%%%%%%%%%%%%%
        
        %%% Distance between centroid %%%
        DistCentroid_MIP = sqrt( (Cent_Im_MIP.Centroid(1) - Cent_Im_Iso.Centroid(1))^2 + (Cent_Im_MIP.Centroid(2) - Cent_Im_Iso.Centroid(2))^2);
        
        DistCentroid_SIP = sqrt( (Cent_Im_SIP.Centroid(1) - Cent_Im_Iso.Centroid(1))^2 + (Cent_Im_SIP.Centroid(2) - Cent_Im_Iso.Centroid(2))^2);
        %%%%%%%%%%%%%%%%%%%%
        
        
%         %%%
%         figure('Name',Animal_Name)
%         subplot(2,3,1)
%         imagesc(Im_MIP_Th)
%         title('MIP thresh')
%         subplot(2,3,2)
%         imagesc(M_Iso_Th)
%         title('Iso thresh')
%         subplot(2,3,3)
%         imagesc(Over_MIP_Iso_Inter)
%         title('Overlapped Area')
%         
%         subplot(2,3,4)
%         imagesc(Im_SIP_Th)
%         title('SIP thresh')
%         subplot(2,3,5)
%         imagesc(M_Iso_Th)
%         title('Iso thresh')
%         subplot(2,3,6)
%         imagesc(Over_SIP_Iso_Inter)
%         title('Overlapped Area')
%         %%%
%         saveas(gca,[Animal_Name,'_',Str_MIP_Week{i_week}],'jpeg');
%         pause
        
        
        
        %%%%%%%%%%%%%%%%%% Store %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strfind(Animal_Name,'BoNT')
            
            if i_week == 1
                %week 1
                ind_TREAT = 3;
            else
                %week 4
                ind_TREAT = 4;
            end
            
        elseif strfind(Animal_Name,'stroke')
            
            ind_TREAT = 2;
        elseif strfind(Animal_Name,'control')
            
            ind_TREAT = 1;
        end
        
        
        Area_Over_MIP_STORE{ ind_TREAT,1}    = [Area_Over_MIP_STORE{ ind_TREAT,1}    Area_Over_MIP];
        Area_Over_SIP_STORE{ ind_TREAT,1}    = [Area_Over_SIP_STORE{ ind_TREAT,1}    Area_Over_SIP];
        
        DistCentroid_MIP_STORE{ ind_TREAT,1} = [DistCentroid_MIP_STORE{ ind_TREAT,1} DistCentroid_MIP];
        DistCentroid_SIP_STORE{ ind_TREAT,1} = [DistCentroid_SIP_STORE{ ind_TREAT,1} DistCentroid_SIP];
        
        
        
        
    end %end for week
    %%%%%%%%%%
    
    
end % LaTO

%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot AREA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MIP_TOT_Mean = [];
MIP_TOT_Std  = [];
SIP_TOT_Mean = [];
SIP_TOT_Std  = [];

for ind_TREAT =1:size(Area_Over_MIP_STORE,1)
    
    for i_MS = 1:2
        
        if i_MS == 1
            %%% MIP %%
            CountData = Area_Over_MIP_STORE{ ind_TREAT,1};
        else
            %%% SIP %%
            CountData = Area_Over_SIP_STORE{ ind_TREAT,1};
        end
        
        %mean animals
%         CountDataNorm_MEAN = nanmedian(CountData,2);
        CountDataNorm_MEAN = nanmean(CountData,2);
        %std animals
        CountDataNorm_STD  = nanstd(CountData,[],2)/sqrt(length(CountData));
        
        if i_MS == 1
            %%% MIP %%
            MIP_TOT_Mean = [MIP_TOT_Mean     CountDataNorm_MEAN];
            MIP_TOT_Std  = [MIP_TOT_Std      CountDataNorm_STD];
        else
            %%% SIP %%
            SIP_TOT_Mean = [SIP_TOT_Mean     CountDataNorm_MEAN];
            SIP_TOT_Std  = [SIP_TOT_Std      CountDataNorm_STD];
        end
        
    end
    
end

figure('Name','Plot')
subplot(1,2,1)
errorbar(1:size(Area_Over_MIP_STORE,1),MIP_TOT_Mean,MIP_TOT_Std)
set(gca,'xTick',     1:size(Area_Over_MIP_STORE,1))
set(gca,'xTickLabel',{'Control', 'Stroke', 'RehabW1','RehabW4'})
rotateXLabels( gca, 60)
title('MIP')
ylabel('overlapped area bw MIP and Iso [px]')

subplot(1,2,2)
errorbar(1:size(Area_Over_SIP_STORE,1),SIP_TOT_Mean, SIP_TOT_Std)
set(gca,'xTick',     1:size(Area_Over_SIP_STORE,1))
set(gca,'xTickLabel',{'Control', 'Stroke', 'RehabW1','RehabW4'})
rotateXLabels( gca, 60)
title('SIP')
ylabel('overlapped area bw SIP and Iso [px]')

filename_1 = ['Overlap_Area_MIPSIP_and_Isochr_FIG'];
% saveas(gca,[filename_1],'fig');
% saveas(gca,[filename_1],'jpeg');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot DISTANCE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MIP_TOT_Mean = [];
MIP_TOT_Std  = [];
SIP_TOT_Mean = [];
SIP_TOT_Std  = [];

for ind_TREAT =1:size(Area_Over_MIP_STORE,1)
    
    for i_MS = 1:2
        
        if i_MS == 1
            %%% MIP %%
            CountData = DistCentroid_MIP_STORE{ ind_TREAT,1};
        else
            %%% SIP %%
            CountData = DistCentroid_SIP_STORE{ ind_TREAT,1};
        end
        
        %mean animals
        CountDataNorm_MEAN = nanmean(CountData,2);
%         CountDataNorm_MEAN = nanmedian(CountData,2);
        %std animals
        CountDataNorm_STD  = nanstd(CountData,[],2);
        
        if i_MS == 1
            %%% MIP %%
            MIP_TOT_Mean = [MIP_TOT_Mean     CountDataNorm_MEAN];
            MIP_TOT_Std  = [MIP_TOT_Std      CountDataNorm_STD];
        else
            %%% SIP %%
            SIP_TOT_Mean = [SIP_TOT_Mean     CountDataNorm_MEAN];
            SIP_TOT_Std  = [SIP_TOT_Std      CountDataNorm_STD];
        end
        
    end
    
end

figure('Name','Plot')
subplot(1,2,1)
errorbar(1:size(Area_Over_MIP_STORE,1),MIP_TOT_Mean,MIP_TOT_Std)
set(gca,'xTick',     1:size(Area_Over_MIP_STORE,1))
set(gca,'xTickLabel',{'Control', 'Stroke', 'RehabW1','RehabW4'})
rotateXLabels( gca, 60)
title('MIP')
ylabel('distance MIP and Iso centroids [px]')

subplot(1,2,2)
errorbar(1:size(Area_Over_SIP_STORE,1),SIP_TOT_Mean,SIP_TOT_Std)
set(gca,'xTick',     1:size(Area_Over_SIP_STORE,1))
set(gca,'xTickLabel',{'Control', 'Stroke', 'RehabW1','RehabW4'})
rotateXLabels( gca, 60)
title('SIP')
ylabel('distance SIP and Iso centroids [px]')

filename_1 = ['CentroidDistance_MIPSIP_vs_Isochr_FIG'];
% saveas(gca,[filename_1],'fig');
% saveas(gca,[filename_1],'jpeg');
