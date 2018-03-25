
User = getenv('USERNAME');

folder = ['C:\Users\',User,'\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\05_Maximum_Intensity_Projection\_Store_MIP_all_Days_and_MAT'];


ListAnimalTogether = {  'GCaMPChR2_7_control','GCaMPChR2_17_control','GCaMPChR2_23_control','GCaMPChR2_24_control',...
                        ...
                        'GCaMPChR2_8_stroke', 'GCaMPChR2_9_stroke', 'GCaMPChR2_19_stroke', 'GCaMPChR2_22_stroke','GCaMPChR2_25_stroke',...
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



for i=1:length(ListAnimalTogether)
    
    
    folderAn = ListAnimalTogether{i};
    

        if strfind(folderAn,'stroke_BoNT')
            %rehab 1
            fileA = load( [folder,'\',folderAn,'\','Im_AlongDays_Week_1_',folderAn,'_MIP_max']);
            MIP_max_REHAB1  = MIP_max_REHAB1+fileA.Im_med/max(max(fileA.Im_med));
            clear fileA
            %rehab 4
            fileA = load( [folder,'\',folderAn,'\','Im_AlongDays_Week_4_',folderAn,'_MIP_max']);
            MIP_max_REHAB4  = MIP_max_REHAB4+fileA.Im_med/max(max(fileA.Im_med));
            clear fileA
            
           
            %rehab 1
            fileA = load( [folder,'\',folderAn,'\','Im_AlongDays_Week_1_',folderAn,'_MIP_sum']);
            MIP_sum_REHAB1  = MIP_sum_REHAB1+fileA.Im_med/max(max(fileA.Im_med));
            clear fileA
            %rehab 4
            fileA = load( [folder,'\',folderAn,'\','Im_AlongDays_Week_4_',folderAn,'_MIP_sum']);
            MIP_sum_REHAB4  = MIP_sum_REHAB4+fileA.Im_med/max(max(fileA.Im_med));
            clear fileA
             
            iter(1) = iter(1)+1;
            iter(2) = iter(2)+1;
            
        elseif strfind(folderAn,'stroke')
            %stroke
            fileA = load( [folder,'\',folderAn,'\','Im_AlongDays_',folderAn,'_MIP_max']);
            MIP_max_STROKE = MIP_max_STROKE+fileA.Im_med/max(max(fileA.Im_med));
            clear fileA
            
            fileA = load( [folder,'\',folderAn,'\','Im_AlongDays_',folderAn,'_MIP_sum']);
            MIP_sum_STROKE  = MIP_sum_STROKE+fileA.Im_med/max(max(fileA.Im_med));
            clear fileA
            
            iter(3) = iter(3)+1;
            
        elseif strfind(folderAn,'control')
            %control
            fileA = load( [folder,'\',folderAn,'\','Im_AlongDays_',folderAn,'_MIP_max']);
            MIP_max_CONTROL = MIP_max_CONTROL+fileA.Im_med/max(max(fileA.Im_med));
            clear fileA
            
            fileA = load( [folder,'\',folderAn,'\','Im_AlongDays_',folderAn,'_MIP_sum']);
            MIP_sum_CONTROL  = MIP_sum_CONTROL+fileA.Im_med/max(max(fileA.Im_med));
            clear fileA
            
            iter(4) = iter(4)+1;
            
        end
    
end

sz = size(MIP_max_REHAB1);

figure('Name','MIP')
subplot(221)
% MIP_max_REHAB1(MIP_max_REHAB1>20) = 20;
imagesc(MIP_max_REHAB1/iter(1))
caxis([0 1])
title('Rehab 1st week')
changeLabel_MIP(sz(1),sz(2),1);
subplot(222)
% MIP_max_REHAB4(MIP_max_REHAB4>20) = 20;
imagesc(MIP_max_REHAB4/iter(2))
caxis([0 1])
title('Rehab 4th week')
changeLabel_MIP(sz(1),sz(2),1);
subplot(223)
% MIP_max_STROKE(MIP_max_STROKE>20) = 20;
imagesc(MIP_max_STROKE/iter(3))
caxis([0 1])
title('Stroke 4th week')
changeLabel_MIP(sz(1),sz(2),1);
subplot(224)
% MIP_max_CONTROL(MIP_max_CONTROL>20) = 20;
imagesc(MIP_max_CONTROL/iter(4))
caxis([0 1])
title('Control')
changeLabel_MIP(sz(1),sz(2),1);
colormap hot

figure('Name','SIP')
subplot(221)
imagesc(MIP_sum_REHAB1/iter(1))
caxis([0 1])
title('Rehab 1st week')
changeLabel_MIP(sz(1),sz(2),1);
subplot(222)
imagesc(MIP_sum_REHAB4/iter(2))
caxis([0 1])
title('Rehab 4th week')
changeLabel_MIP(sz(1),sz(2),1);
subplot(223)
imagesc(MIP_sum_STROKE/iter(3))
caxis([0 1])
title('Stroke 4th week')
changeLabel_MIP(sz(1),sz(2),1);
changeLabel_MIP(sz(1),sz(2),1);
subplot(224)
imagesc(MIP_sum_CONTROL/iter(4))
caxis([0 1])
title('Control')
changeLabel_MIP(sz(1),sz(2),1);
colormap hot






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
