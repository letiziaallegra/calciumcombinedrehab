% -carico SIP della settimana (valore 0-5)
% -divido per il valore max (in modo tale da evitare casi in cui ho meno di 5
%     giorni per settimana)
% -vado a vedere dove c'è attivazione per almeno 3/5 giorni
% -queste SIP sogliate le vado a mettere tutte insieme e vedo in media a vedere cosa
%  succede -> uso un'ultreiore soglia di 3/5

close all
clear all
clc

User = getenv('USERNAME');

folder = ['C:\Users\',User,'\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\05_Maximum_Intensity_Projection_SeqSelez\_Store_MIP_SIP_all_Days_and_MAT'];


% ListAnimalTogether = {  'GCaMPChR2_7_control','GCaMPChR2_17_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
%                         ...
%                         'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke','GCaMPChR2_25_stroke','GCaMPChR2_26_stroke'};

                    ListAnimalTogether = {      'GCaMPChR2_7_control','GCaMPChR2_17_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
                                                'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke','GCaMPChR2_25_stroke', 'GCaMPChR2_26_stroke'...
                                                'GCaMP16_stroke_BoNT','GCaMP18_stroke_BoNT',...
                                                'GCaMPChR2_11_stroke_BoNT', 'GCaMPChR2_12_stroke_BoNT',...
                                                'GCaMPChR2_14_stroke_BoNT', 'GCaMPChR2_15_stroke_BoNT'};
                    
                    
MIP_max_REHAB1 = zeros(512,512);
MIP_max_REHAB4 = zeros(512,512);
MIP_max_STROKE  = zeros(512,512);
MIP_max_CONTROL = zeros(512,512);

MIP_sum_REHAB1 = zeros(512,512);
MIP_sum_REHAB4 = zeros(512,512);
MIP_sum_STROKE  = zeros(512,512);
MIP_sum_CONTROL = zeros(512,512);

iter = zeros(4,1);

Th_single = 3/5;

for i=1:length(ListAnimalTogether)
    
    
    folderAn = ListAnimalTogether{i};
    

        if strfind(folderAn,'stroke_BoNT')
            %rehab 1
            fileA = load( [folder,'\',folderAn,'\','Im_AlongDays_Week_1_',folderAn,'_Long_SelSeq_MIP']);
            MIP_max_REHAB1  = MIP_max_REHAB1+fileA.Im_med/max(max(fileA.Im_med));
            clear fileA
            %rehab 4
            fileA = load( [folder,'\',folderAn,'\','Im_AlongDays_Week_4_',folderAn,'_Long_SelSeq_MIP']);
            MIP_max_REHAB4  = MIP_max_REHAB4+fileA.Im_med/max(max(fileA.Im_med));
            clear fileA
            
           
            %rehab 1
            fileA = load( [folder,'\',folderAn,'\','Im_AlongDays_Week_1_',folderAn,'_Long_SelSeq_SIP']);
            imagesc(double(fileA.Im_med/max(max(fileA.Im_med))>Th_single))
            MIP_sum_REHAB1  = MIP_sum_REHAB1+double(fileA.Im_med/max(max(fileA.Im_med))>Th_single);
            clear fileA
            %rehab 4
            fileA = load( [folder,'\',folderAn,'\','Im_AlongDays_Week_4_',folderAn,'_Long_SelSeq_SIP']);
            imagesc(double(fileA.Im_med/max(max(fileA.Im_med))>Th_single))
            MIP_sum_REHAB4  = MIP_sum_REHAB4+double(fileA.Im_med/max(max(fileA.Im_med))>Th_single);
            clear fileA
             
            iter(1) = iter(1)+1;
            iter(2) = iter(2)+1;
            
        elseif strfind(folderAn,'stroke')
            %stroke
            fileA = load( [folder,'\',folderAn,'\','Im_AlongDays_',folderAn,'_Long_SelSeq_MIP']);
            MIP_max_STROKE = MIP_max_STROKE+fileA.Im_med/max(max(fileA.Im_med));
            clear fileA
            
            fileA = load( [folder,'\',folderAn,'\','Im_AlongDays_',folderAn,'_Long_SelSeq_SIP']);
            imagesc(double(fileA.Im_med/max(max(fileA.Im_med))>Th_single))
            MIP_sum_STROKE  = MIP_sum_STROKE+double(fileA.Im_med/max(max(fileA.Im_med))>Th_single);
            clear fileA
            
            iter(3) = iter(3)+1;
            
        elseif strfind(folderAn,'control')
            %control
            fileA = load( [folder,'\',folderAn,'\','Im_AlongDays_',folderAn,'_Long_SelSeq_MIP']);
            MIP_max_CONTROL = MIP_max_CONTROL+fileA.Im_med/max(max(fileA.Im_med));
            clear fileA
            
            fileA = load( [folder,'\',folderAn,'\','Im_AlongDays_',folderAn,'_Long_SelSeq_SIP']);
            imagesc(double(fileA.Im_med/max(max(fileA.Im_med))>Th_single))
            MIP_sum_CONTROL  = MIP_sum_CONTROL+double(fileA.Im_med/max(max(fileA.Im_med))>Th_single);
            clear fileA
            
            iter(4) = iter(4)+1;
            
        end
    
end

sz = size(MIP_max_REHAB1);

% figure('Name','MIP')
% subplot(223)
% % MIP_max_REHAB1(MIP_max_REHAB1>20) = 20;
% imagesc(MIP_max_REHAB1/iter(1))
% caxis([0 1])
% title('Rehab 1st week')
% changeLabel_MIP(sz(1),sz(2),1);
% 
% subplot(224)
% % MIP_max_REHAB4(MIP_max_REHAB4>20) = 20;
% imagesc(MIP_max_REHAB4/iter(2))
% caxis([0 1])
% title('Rehab 4th week')
% changeLabel_MIP(sz(1),sz(2),1);
% 
% subplot(222)
% % MIP_max_STROKE(MIP_max_STROKE>20) = 20;
% imagesc(MIP_max_STROKE/iter(3))
% caxis([0 1])
% title('Stroke 4th week')
% changeLabel_MIP(sz(1),sz(2),1);
% 
% subplot(221)
% % MIP_max_CONTROL(MIP_max_CONTROL>20) = 20;
% imagesc(MIP_max_CONTROL/iter(4))
% caxis([0 1])
% title('Control')
% changeLabel_MIP(sz(1),sz(2),1);
% colormap hot


%%% -> N.B. HO SOMMATO LE SIP e poi divise per il loro numero (quindi ho fatto la media)
figure('Name','SIP') 
subplot(223)
imagesc(MIP_sum_REHAB1/iter(1))
caxis([0 1])
title('Rehab 1st week')
changeLabel_MIP(sz(1),sz(2),1);
subplot(224)
imagesc(MIP_sum_REHAB4/iter(2))
caxis([0 1])
title('Rehab 4th week')
changeLabel_MIP(sz(1),sz(2),1);
subplot(222)
imagesc(MIP_sum_STROKE/iter(3))
caxis([0 1])
title('Stroke 4th week')
changeLabel_MIP(sz(1),sz(2),1);
changeLabel_MIP(sz(1),sz(2),1);
subplot(221)
imagesc(MIP_sum_CONTROL/iter(4))
caxis([0 1])
title('Control')
changeLabel_MIP(sz(1),sz(2),1);
colormap hot


%%% THRESHOLDED %%% -> N.B. HO SOMMATO LE SIP e poi divise per il loro numero (quindi ho fatto la media)
THS = 3/5;
figSIP_T = figure('Name','SIP')
subplot(223)
MIP_sum_REHAB1_T = MIP_sum_REHAB1/iter(1);
MIP_sum_REHAB1_T(MIP_sum_REHAB1_T<=THS) = 0;
imagesc(MIP_sum_REHAB1_T)
caxis([0 1])
title('Rehab 1st week')
changeLabel_MIP(sz(1),sz(2),1);

subplot(224)
MIP_sum_REHAB4_T = MIP_sum_REHAB4/iter(2);
MIP_sum_REHAB4_T(MIP_sum_REHAB4_T<=THS) = 0;
imagesc(MIP_sum_REHAB4_T)
caxis([0 1])
title('Rehab 4th week')
changeLabel_MIP(sz(1),sz(2),1);

subplot(222)
MIP_sum_STROKE_T = MIP_sum_STROKE/iter(3);
MIP_sum_STROKE_T(MIP_sum_STROKE_T<=THS) = 0;
imagesc(MIP_sum_STROKE_T)
caxis([0 1])
title('Stroke 4th week')
changeLabel_MIP(sz(1),sz(2),1);

subplot(221)
MIP_sum_CONTROL_T = MIP_sum_CONTROL/iter(3);
MIP_sum_CONTROL_T(MIP_sum_CONTROL_T<=THS) = 0;
imagesc(MIP_sum_CONTROL_T)
caxis([0 1])
title('Control')
changeLabel_MIP(sz(1),sz(2),1);
colormap hot
%%%%%%%%%%%%%%%%%%



%%%% calcola il centro di massa e l'area sogliata per il caso precedente %%
%AREA
RappPixelmm = (4.4/512); % from px to mm
Area_Pixel = RappPixelmm^2;% [mm2] -> quadrato totale di 4.4mm di lato disposti su 512 px di lato
MIP_sum_CONTROL_T_1 = MIP_sum_CONTROL_T>0;
MIP_sum_STROKE_T_1  = MIP_sum_STROKE_T>0;
MIP_sum_REHAB1_T_1  = MIP_sum_REHAB1_T>0;
MIP_sum_REHAB4_T_1  = MIP_sum_REHAB4_T>0;

Area_CNT = sum(sum(double(MIP_sum_CONTROL_T_1)))*Area_Pixel;
Area_STR = sum(sum(double(MIP_sum_STROKE_T_1)))*Area_Pixel;
Area_R1W = sum(sum(double(MIP_sum_REHAB1_T_1)))*Area_Pixel;
Area_R4W = sum(sum(double(MIP_sum_REHAB4_T_1)))*Area_Pixel;
Area = [Area_CNT, Area_STR, Area_R1W, Area_R4W];

%Centroide
[y x] = find( MIP_sum_CONTROL_T_1 );
cent = [mean(x) mean(y)];
dCenCNT        = sqrt(cent(1)^2+cent(2)^2)* RappPixelmm;
figure(figSIP_T); subplot(221); hold on; scatter(cent(1), cent(2),'r')

[y x] = find( MIP_sum_STROKE_T_1 );
cent = [mean(x) mean(y)];
dCenSTR         = sqrt(cent(1)^2+cent(2)^2)* RappPixelmm;
figure(figSIP_T); subplot(222); hold on; scatter(cent(1), cent(2),'r')

[y x] = find( MIP_sum_REHAB1_T_1 );
cent = [mean(x) mean(y)];
dCenR1W        = sqrt(cent(1)^2+cent(2)^2)* RappPixelmm;
figure(figSIP_T); subplot(223); hold on; scatter(cent(1), cent(2),'r')

[y x] = find( MIP_sum_REHAB4_T_1 );
cent = [mean(x) mean(y)];
dCenR4W         = sqrt(cent(1)^2+cent(2)^2)* RappPixelmm;
figure(figSIP_T); subplot(224); hold on; scatter(cent(1), cent(2),'r')

dCentroid = [dCenCNT, dCenSTR, dCenR1W, dCenR4W];

figure
scatter(dCentroid,Area)
text(dCentroid+0.05, Area+0.05,{'Control','Stroke','Rehab1w','Rehab4w'})
xlabel('Distance from Bregma [mm]')
ylabel('Area [mm_{^2}]')




% figure('Name','MIP')
% subplot(221)
% % MIP_max_REHAB1(MIP_max_REHAB1>20) = 20;
% imshow(MIP_max_REHAB1/iter(1))
% caxis([0 1])
% title('Rehab 1st week')
% changeLabel_MIP(sz(1),sz(2),1);
% subplot(222)
% % MIP_max_REHAB4(MIP_max_REHAB4>20) = 20;
% imshow(MIP_max_REHAB4/iter(2))
% caxis([0 1])
% title('Rehab 4th week')
% changeLabel_MIP(sz(1),sz(2),1);
% subplot(223)
% % MIP_max_STROKE(MIP_max_STROKE>20) = 20;
% imshow(MIP_max_STROKE/iter(3))
% caxis([0 1])
% title('Stroke 4th week')
% changeLabel_MIP(sz(1),sz(2),1);
% subplot(224)
% % MIP_max_CONTROL(MIP_max_CONTROL>20) = 20;
% imshow(MIP_max_CONTROL/iter(4))
% caxis([0 1])
% title('Control')
% changeLabel_MIP(sz(1),sz(2),1);
% colormap hot
% 
% figure('Name','SIP')
% subplot(221)
% imshow(MIP_sum_REHAB1/iter(1))
% caxis([0 1])
% title('Rehab 1st week')
% changeLabel_MIP(sz(1),sz(2),1);
% subplot(222)
% imshow(MIP_sum_REHAB4/iter(2))
% caxis([0 1])
% title('Rehab 4th week')
% changeLabel_MIP(sz(1),sz(2),1);
% subplot(223)
% imshow(MIP_sum_STROKE/iter(3))
% caxis([0 1])
% title('Stroke 4th week')
% changeLabel_MIP(sz(1),sz(2),1);
% changeLabel_MIP(sz(1),sz(2),1);
% subplot(224)
% imshow(MIP_sum_CONTROL/iter(4))
% caxis([0 1])
% title('Control')
% changeLabel_MIP(sz(1),sz(2),1);
% colormap hot
% 
% 
% 
% 
