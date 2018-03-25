%
% Find percentage of M1 and other areas that constitute the ROI of the MIP
%


clear
close all
clc

CurrDir = cd;


%%%% Choice of the animal and trial day
UsbPortHD = 'I';
UserName  = getenv('username');
AnimalDir = [UsbPortHD,':\LENS\Animals Data'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Animal_Name %%%%%%%%%%%%%%%%%%%
% ListAnimalTogether = { 'GCaMPChR2_7_control'};


ListAnimalTogether = {  'GCaMPChR2_7_control','GCaMPChR2_24_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
                        'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke','GCaMPChR2_25_stroke', 'GCaMPChR2_26_stroke',...
                        'GCaMP16_stroke_BoNT','GCaMP18_stroke_BoNT',...
                        'GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT',...
                        'GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT'};
                 
% ListAnimalTogether = {  'GCaMPChR2_7_control','GCaMPChR2_24_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
%                         'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke','GCaMPChR2_25_stroke',...
%                         ...
%                         'GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT',...
%                         'GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ThreshMIP = 3;


%%%%%%%% load the list of indeces corresponding to the functional areas%
load('RegionArea_FunctionalMask_Index'); %-> RegionArea_FunctionalMask_Index
AREA_NAME   = {'M2','CFA','M1Sh','M1HL','M1Tk','S1','S1Bc','RSD','LPTa','V1','AuD'};
FuncAres_ToTake      = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
% AREA_NAME   = {'CFA','M1Sh','M1HL','M1Tk','S1','S1Bc','RSD','LPTa','V1'};
% FuncAres_ToTake      = {2, 3, 4, 5, 6, 7, 8, 9, 10};
% AREA_NAME   = {'M1','S1','S1BF','RSD','V1'};
% FuncAres_ToTake      = {  2, [3 4 6], 5, 9, 10};
Num_FunAreas  = size(RegionArea_FunctionalMask_Index,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% list animal for plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_N_control   = 'GCaMPChR2_7_control_Week_1';
P_N_stroke    = 'GCaMPChR2_25_stroke_Week_1';
P_N_rehab_1W  = 'GCaMPChR2_14_stroke_BoNT_Week_1';
P_N_rehab_4W  = 'GCaMPChR2_14_stroke_BoNT_Week_4';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%count index in M ROI
Count_Mk_MIP_STORE = cell(4,1);
Count_Mk_SIP_STORE = cell(4,1);


for LATo = 1:length(ListAnimalTogether)
    
%     clearvars -except CurrDir UsbPortHD UserName AnimalDir ListAnimalTogether LATo
    
    Animal_Name                    = [ListAnimalTogether{LATo}];
    
    %%%%%%%% TrialDay %%%%%%%%%%%%%%%%%%%%%%%
    % TrialDay_choice = '01';
    TrialDay_choice = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% LOAD data REGION BOUND MIP %%%%%%%%%%%%%%%%%%%%
    ROI_Folder   = ['C:\Users\',UserName,'\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\05_Maximum_Intensity_Projection_SeqSelez\_Store_MIP_SIP_all_Days_and_MAT\',Animal_Name];
    
    %%%%%%%%%%
    if strfind(Animal_Name,'BoNT')
        
        %week1
        load([ROI_Folder,'\','Im_AlongDays_Week_1_',Animal_Name,'_Long_SelSeq_MIP']);
        Im_MIP = Im_med;
        Im_MIP = Im_MIP>=ThreshMIP;        
        clear Im_med
        load([ROI_Folder,'\','Im_AlongDays_Week_1_',Animal_Name,'_Long_SelSeq_SIP']);
        Im_SIP = Im_med;
        Im_SIP = Im_SIP>=ThreshMIP;
        clear Im_med
                
        %week4
        load([ROI_Folder,'\','Im_AlongDays_Week_4_',Animal_Name,'_Long_SelSeq_MIP']);
        Im_MIP_4 = Im_med;
        Im_MIP_4 = Im_MIP_4>=ThreshMIP; 
        clear Im_med
        load([ROI_Folder,'\','Im_AlongDays_Week_4_',Animal_Name,'_Long_SelSeq_SIP']);
        Im_SIP_4 = Im_med;
        Im_SIP_4 = Im_SIP_4>=ThreshMIP; 
        clear Im_med
        
        Check_Week = 2;
        
    else
        
        load([ROI_Folder,'\','Im_AlongDays_',Animal_Name,'_Long_SelSeq_MIP']);
        Im_MIP = Im_med;
        Im_MIP = Im_MIP>=ThreshMIP;        
        clear Im_med
        load([ROI_Folder,'\','Im_AlongDays_',Animal_Name,'_Long_SelSeq_SIP']);
        Im_SIP = Im_med;
        Im_SIP = Im_SIP>=ThreshMIP;
        clear Im_med
        
        Check_Week = 1;
        
    end
    %%%%%%%%%%          
   
    [s1 s2] = size(Im_MIP);
    Count_Mk_MIP = [];
    Count_Mk_SIP = [];
    Count_Mk_MIP_4 = [];
    Count_Mk_SIP_4 = [];
    
    for i_A=1:length(FuncAres_ToTake)
        
        for i_ch=1:Check_Week %for check
            
            Area_index = FuncAres_ToTake{i_A};
            
            MaskArea = zeros(s1,s2);
            MaskArea(sub2ind([s1 s2],RegionArea_FunctionalMask_Index{Area_index ,1}(:,2),RegionArea_FunctionalMask_Index{Area_index ,1}(:,1))) = 1;
            
            if i_ch==1 %if check
                
                %Part of MIP area and SIP area that are in the M area
                Mk_MIP = MaskArea.*Im_MIP;
                Mk_SIP = MaskArea.*Im_SIP;
                
                %count number of points in the M area
                Count_Mk_MIP = [Count_Mk_MIP; sum(sum(Mk_MIP))];
                Count_Mk_SIP = [Count_Mk_SIP; sum(sum(Mk_SIP))];
                
%                 %%%
%                 figure
%                 subplot(2,3,[1 4])
%                 imagesc(MaskArea)
%                 subplot(2,3,2)
%                 imagesc(Im_MIP)
%                 subplot(2,3,5)
%                 imagesc(Mk_MIP)
%                 subplot(2,3,3)
%                 imagesc(Im_SIP)
%                 subplot(2,3,6)
%                 imagesc(Mk_SIP)
%                 ciao = 1;
%                 %%%
                
            else
                
                %Part of MIP area and SIP area that are in the M area
                Mk_MIP_4 = MaskArea.*Im_MIP_4;
                Mk_SIP_4 = MaskArea.*Im_SIP_4;
                
                %count number of points in the M area
                Count_Mk_MIP_4 = [Count_Mk_MIP_4; sum(sum(Mk_MIP_4))];
                Count_Mk_SIP_4 = [Count_Mk_SIP_4; sum(sum(Mk_SIP_4))];
                
%                 %%%
%                 figure
%                 subplot(2,3,[1 4])
%                 imagesc(MaskArea)
%                 subplot(2,3,2)
%                 imagesc(Im_MIP_4)
%                 subplot(2,3,5)
%                 imagesc(Mk_MIP_4)
%                 subplot(2,3,3)
%                 imagesc(Im_SIP_4)
%                 subplot(2,3,6)
%                 imagesc(Mk_SIP_4)
%                 ciao = 1;
%                 %%%
                
            end %end if check
            
        end %end for check
        
    end
    
    
    
    if strfind(Animal_Name,'BoNT')
        %week 4
        ind_TREAT = 4;
        Count_Mk_MIP_STORE{ ind_TREAT,1} = [Count_Mk_MIP_STORE{ ind_TREAT,1} Count_Mk_MIP_4];
        Count_Mk_SIP_STORE{ ind_TREAT,1} = [Count_Mk_SIP_STORE{ ind_TREAT,1} Count_Mk_SIP_4];
        %week 1
        ind_TREAT = 3;
        
    elseif strfind(Animal_Name,'stroke')        
        ind_TREAT = 2;
    elseif strfind(Animal_Name,'control')       
        ind_TREAT = 1;
    end
    
    Count_Mk_MIP_STORE{ ind_TREAT,1} = [Count_Mk_MIP_STORE{ ind_TREAT,1} Count_Mk_MIP];
    Count_Mk_SIP_STORE{ ind_TREAT,1} = [Count_Mk_SIP_STORE{ ind_TREAT,1} Count_Mk_SIP];
    
        
end % LaTO 
    
%



%%% plot %%%
CountDataNorm_MEAN_TOT_MIP = [];
CountDataNorm_MEAN_TOT_SIP = [];
for ind_TREAT =1:size(Count_Mk_MIP_STORE,1)
    
    for i_MS = 1:2
        
        if i_MS == 1
            %%% MIP %%
            CountData = Count_Mk_MIP_STORE{ ind_TREAT,1};
        else
            %%% SIP %%
            CountData = Count_Mk_SIP_STORE{ ind_TREAT,1};
        end

        %sum single animal
        SumData   = sum(CountData,1);
        
        %percentage single animal 
        CountDataNorm = CountData ./ repmat(SumData,size(CountData,1),1) *100;
        
        %mean animals
        CountDataNorm_MEAN = nanmean(CountDataNorm,2);
        %std animals
        CountDataNorm_STD  = nanstd(CountDataNorm,[],2);

        if i_MS == 1
            %%% MIP %%
            CountDataNorm_MEAN_TOT_MIP = [CountDataNorm_MEAN_TOT_MIP CountDataNorm_MEAN];
        else
            %%% SIP %%
            CountDataNorm_MEAN_TOT_SIP = [CountDataNorm_MEAN_TOT_SIP CountDataNorm_MEAN];
        end
                       
    end
    
end

figure('Name','Plot')
subplot(1,2,1)
bar(1:size(Count_Mk_MIP_STORE,1),CountDataNorm_MEAN_TOT_MIP','stacked')
set(gca,'xTick',     1:size(Count_Mk_MIP_STORE,1))
set(gca,'xTickLabel',{'Control', 'Stroke', 'RehabW1','RehabW4'})
rotateXLabels( gca, 60)
title('MIP')
ylabel('% of areas constituting MIP')
ylim([0 100])

subplot(1,2,2)
bar(1:size(Count_Mk_SIP_STORE,1),CountDataNorm_MEAN_TOT_SIP','stacked')
set(gca,'xTick',     1:size(Count_Mk_SIP_STORE,1))
set(gca,'xTickLabel',{'Control', 'Stroke', 'RehabW1','RehabW4'})
rotateXLabels( gca, 60)
title('SIP')
ylabel('% of areas constituting SIP')
ylim([0 100])

legend(AREA_NAME)
colormap hot

