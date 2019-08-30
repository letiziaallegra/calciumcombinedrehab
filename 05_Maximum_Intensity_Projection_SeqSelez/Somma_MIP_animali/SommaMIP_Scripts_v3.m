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

%folder = ['C:\Users\',User,'\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\05_Maximum_Intensity_Projection_SeqSelez\_Store_MIP_SIP_all_Days_and_MAT'];
folder = ['C:\Users\CNR-SSSUP\Google Drive_SL\Piattaforma Stefano\ELABORAZIONE DATA\05_Maximum_Intensity_Projection_SeqSelez\_Store_MIP_SIP_all_Days_and_MAT'];
folder = ['/Users/alessandro/Desktop/ELABORAZIONE DATA/05_Maximum_Intensity_Projection_SeqSelez/_Store_MIP_SIP_all_Days_and_MAT/'];

% ListAnimalTogether = {  'GCaMPChR2_7_control','GCaMPChR2_17_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
%                         ...
%                         'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke','GCaMPChR2_25_stroke','GCaMPChR2_26_stroke'};

% LIST OF FILTERED ANIMALS
ListAnimalTogether = {      'GCaMPChR2_7_control','GCaMPChR2_17_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
    'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke','GCaMPChR2_25_stroke', 'GCaMPChR2_26_stroke'...
    'GCamp_25_sani','GCamp_26_sani'};


MIP_max_REHAB1 = zeros(103,103);
MIP_max_REHAB4 = zeros(103,103);
MIP_max_STROKE  = zeros(103,103);
MIP_max_CONTROL = zeros(103,103);
MIP_max_STIM = zeros(100,100);
MIP_max_SANI = zeros(100,100);


MIP_sum_REHAB1 = zeros(103,103);
MIP_sum_REHAB4 = zeros(103,103);
MIP_sum_STROKE  = zeros(103,103);
MIP_sum_CONTROL = zeros(103,103);
MIP_sum_STIM = zeros(100,100);
MIP_sum_SANI = zeros(100,100);


iter = zeros(5,1);

Th_single = 2/5;
RappPixelmm = (4.4/103); % from px to mm
RappPixelmmSTIM = (5.25/100); % from px to mm
factor = 0.2; %circa 100/512 perchè le immagini STIM sono differenti nel numero di pixel e nelle dimensioni della finestra (100x100 che corrispondono a 5.25x5.25mm)
factordiv=round(100/5.25*4.4); %fattore per modificare le size


%%%%%%%%
%Area
Area_CN        = [];
Area_ST        = [];
Area_R1        = [];
Area_R4        = [];
%Area_S4        = [];
%Centroid
dCen_CN        = [];
dCen_ST        = [];
dCen_R1        = [];
dCen_R4        = [];
%dCen_S4        = [];
%%%%%%%%


for i=1:length(ListAnimalTogether)
    
    
    folderAn = ListAnimalTogether{i};
    
    
    if strfind(folderAn,'stroke_BoNT')
        %rehab 1
        fileA = load( [folder,filesep,folderAn,filesep,'Im_AlongDays_Week_1_',folderAn,'_Long_SelSeq_MIP']);
        MIP_max_REHAB1  = MIP_max_REHAB1+fileA.Im_med/max(max(fileA.Im_med));
        clear fileA
        %rehab 4
        fileA = load( [folder,filesep,folderAn,filesep,'Im_AlongDays_Week_4_',folderAn,'_Long_SelSeq_MIP']);
        MIP_max_REHAB4  = MIP_max_REHAB4+fileA.Im_med/max(max(fileA.Im_med));
        clear fileA
        
        
        %rehab 1
        fileA = load( [folder,filesep,folderAn,filesep,'Im_AlongDays_Week_1_',folderAn,'_Long_SelSeq_SIP']);
        imagesc(double(fileA.Im_med/max(max(fileA.Im_med))>Th_single))
        MIP_sum_REHAB1  = MIP_sum_REHAB1+double(fileA.Im_med/max(max(fileA.Im_med))>Th_single);
        
        %%%%%%%%
        %Area
        fileAT         = (fileA.Im_med/max(max(fileA.Im_med))>Th_single)>0;
        Area_R1        = [Area_R1; sum(sum(double(fileAT)))*RappPixelmm^2];
        %Centroid
        [y x]          = find( double(fileAT) );
        cent           = [mean(x) mean(y)];
        
        dCen_R1        = [dCen_R1;  cent(2)*RappPixelmm];
        %             dCen_R1        = [dCen_R1;sqrt(cent(1)^2+cent(2)^2)*RappPixelmm];
        %%%%%%%%
        clear fileA
        
        
        %rehab 4
        fileA = load( [folder,filesep,folderAn,filesep,'Im_AlongDays_Week_4_',folderAn,'_Long_SelSeq_SIP']);
        imagesc(double(fileA.Im_med/max(max(fileA.Im_med))>Th_single))
        MIP_sum_REHAB4  = MIP_sum_REHAB4+double(fileA.Im_med/max(max(fileA.Im_med))>Th_single);
        
        %%%%%%%%
        %Area
        fileAT         = (fileA.Im_med/max(max(fileA.Im_med))>Th_single)>0;
        Area_R4        = [Area_R4; sum(sum(double(fileAT)))*RappPixelmm^2];
        %Centroid
        [y x]          = find( double(fileAT) );
        cent           = [mean(x) mean(y)];
        
        dCen_R4        = [dCen_R4;  cent(2)*RappPixelmm];
        %             dCen_R4        = [dCen_R4;sqrt(cent(1)^2+cent(2)^2)*RappPixelmm];
        %%%%%%%%
        clear fileA
        
        iter(1) = iter(1)+1;
        iter(2) = iter(2)+1;
        
    elseif strfind(folderAn,'stroke')
        %stroke
        fileA = load( [folder,filesep,folderAn,filesep,'Im_AlongDays_',folderAn,'_Long_SelSeq_MIP']);
        MIP_max_STROKE = MIP_max_STROKE+fileA.Im_med/max(max(fileA.Im_med));
        clear fileA
        
        fileA = load( [folder,filesep,folderAn,filesep,'Im_AlongDays_',folderAn,'_Long_SelSeq_SIP']);
        imagesc(double(fileA.Im_med/max(max(fileA.Im_med))>Th_single))
        MIP_sum_STROKE  = MIP_sum_STROKE+double(fileA.Im_med/max(max(fileA.Im_med))>Th_single);
        
        %%%%%%%%
        %Area
        fileAT         = (fileA.Im_med/max(max(fileA.Im_med))>Th_single)>0;
        Area_ST        = [Area_ST; sum(sum(double(fileAT)))*RappPixelmm^2];
        %Centroid
        [y x]          = find( double(fileAT) );
        cent           = [mean(x) mean(y)];
        
        dCen_ST        = [dCen_ST;  cent(2)*RappPixelmm];
        %             dCen_ST        = [dCen_ST;sqrt(cent(1)^2+cent(2)^2)*RappPixelmm];
        %%%%%%%%
        clear fileA
        
        iter(3) = iter(3)+1;
        
    elseif strfind(folderAn,'control')
        %control
        fileA = load( [folder,filesep,folderAn,filesep,'Im_AlongDays_',folderAn,'_Long_SelSeq_MIP']);
        MIP_max_CONTROL = MIP_max_CONTROL+fileA.Im_med/max(max(fileA.Im_med));
        clear fileA
        
        fileA = load( [folder,filesep,folderAn,filesep,'Im_AlongDays_',folderAn,'_Long_SelSeq_SIP']);
        imagesc(double(fileA.Im_med/max(max(fileA.Im_med))>Th_single))
        MIP_sum_CONTROL  = MIP_sum_CONTROL+double(fileA.Im_med/max(max(fileA.Im_med))>Th_single);
        
        %%%%%%%%
        %Area
        fileAT         = (fileA.Im_med/max(max(fileA.Im_med))>Th_single)>0;
        Area_CN        = [Area_CN; sum(sum(double(fileAT)))*RappPixelmm^2];
        %Centroid
        [y x]          = find( double(fileAT) );
        cent           = [mean(x) mean(y)];
        dCen_CN        = [dCen_CN;  cent(2)*RappPixelmm];
        %             dCen_CN        = [dCen_CN;sqrt(cent(1)^2+cent(2)^2)*RappPixelmm];
        %%%%%%%%
        clear fileA
        
        iter(4) = iter(4)+1;
        
    elseif strfind(folderAn,'sani')
        factor=1;
        fileA = load( [folder,filesep,folderAn,filesep,'Im_AlongDays_',folderAn,'_Long_SelSeq_MIP']);
        MIP_resize=imresize(fileA.Im_med/max(max(fileA.Im_med)),factor);
        MIP_max_SANI = MIP_max_SANI+MIP_resize;
        clear fileA
        
        fileA = load( [folder,filesep,folderAn,filesep,'Im_AlongDays_',folderAn,'_Long_SelSeq_SIP']);
        MIP_sum_resize=imresize(fileA.Im_med/max(max(fileA.Im_med)),factor);
        imagesc(double(MIP_sum_resize>Th_single))
        MIP_sum_SANI  = MIP_sum_SANI+double(MIP_sum_resize>Th_single);
        
        %%%%%%%%
        %Area
        fileAT         = (MIP_sum_resize>Th_single)>0;
        Area_ST        = [Area_ST; sum(sum(double(fileAT)))*RappPixelmmSTIM^2];
        %Centroid
        [y x]          = find( double(fileAT) );
        cent           = [mean(x) mean(y)];
        
        dCen_ST        = [dCen_ST;  cent(2)*RappPixelmmSTIM];
        %             dCen_ST        = [dCen_ST;sqrt(cent(1)^2+cent(2)^2)*RappPixelmmSTIM];
        %%%%%%%%
        clear fileA
        
        iter(2) = iter(2)+1;
        factor = 0.2;
        
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


% %%% -> N.B. HO SOMMATO LE SIP e poi divise per il loro numero (quindi ho fatto la media)
% figure('Name','SIP')
% subplot(223)
% imagesc(MIP_sum_REHAB1/iter(1))
% caxis([0 1])
% title('Rehab 1st week')
% changeLabel_MIP(sz(1),sz(2),1);
% subplot(224)
% imagesc(MIP_sum_REHAB4/iter(2))
% caxis([0 1])
% title('Rehab 4th week')
% changeLabel_MIP(sz(1),sz(2),1);
% subplot(222)
% imagesc(MIP_sum_STROKE/iter(3))
% caxis([0 1])
% title('Stroke 4th week')
% changeLabel_MIP(sz(1),sz(2),1);
% changeLabel_MIP(sz(1),sz(2),1);
% subplot(221)
% imagesc(MIP_sum_CONTROL/iter(4))
% caxis([0 1])
% title('Control')
% changeLabel_MIP(sz(1),sz(2),1);
% colormap hot


%%% THRESHOLDED %%% -> N.B. HO SOMMATO LE SIP e poi divise per il loro numero (quindi ho fatto la media)
% THS = 3/5;
THS = 1/5;
figSIP_T = figure('Name','SIP')
subplot(223)
MIP_sum_REHAB1_T = MIP_sum_REHAB1/iter(1);
MIP_sum_REHAB1_T(MIP_sum_REHAB1_T<=THS) = 0;
imagesc(MIP_sum_REHAB1_T)
caxis([0 1])
% title('Rehab 1st week')
changeLabel_MIP(sz(1),sz(2),1);

subplot(224)
MIP_sum_REHAB4_T = MIP_sum_REHAB4/iter(2);
MIP_sum_REHAB4_T(MIP_sum_REHAB4_T<=THS) = 0;
imagesc(MIP_sum_REHAB4_T)
caxis([0 1])
% title('Rehab 4th week')
changeLabel_MIP(sz(1),sz(2),1);

subplot(222)
MIP_sum_STROKE_T = MIP_sum_STROKE/iter(3);
MIP_sum_STROKE_T(MIP_sum_STROKE_T<=THS) = 0;
imagesc(MIP_sum_STROKE_T)
caxis([0 1])
% title('Stroke 4th week')
changeLabel_MIP(sz(1),sz(2),1);

subplot(221)
MIP_sum_CONTROL_T = MIP_sum_CONTROL/iter(3);
MIP_sum_CONTROL_T(MIP_sum_CONTROL_T<=THS) = 0;
imagesc(MIP_sum_CONTROL_T)
caxis([0 1])
% title('Control')
changeLabel_MIP(sz(1),sz(2),1);

colormap jet

%%%%%%%%%%%%%%%%%%



%%%% calcola il centro di massa e l'area sogliata per il caso precedente %%
%AREA
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


hold on
%Area (dai singoli animali)
mean_Area_SingAn       = [median(Area_CN), median(Area_ST), median(Area_R1), median(Area_R4)];
std_Area_SingAn        = [std(Area_CN), std(Area_ST), std(Area_R1), std(Area_R4)];
mean_dCentroid_SingAn  = [median(dCen_CN), median(dCen_ST), median(dCen_R1), median(dCen_R4)];
std_dCentroid_SingAn   = [std(dCen_CN), std(dCen_ST), std(dCen_R1), std(dCen_R4)];

% scatter(mean_dCentroid_SingAn,mean_Area_SingAn,'r')
figure
errorbarxy(mean_dCentroid_SingAn, mean_Area_SingAn , std_dCentroid_SingAn, std_dCentroid_SingAn, std_Area_SingAn, std_Area_SingAn)
text(mean_dCentroid_SingAn+0.05,mean_Area_SingAn+0.05, {'Control','Stroke','Rehab1w','Rehab4w'})
