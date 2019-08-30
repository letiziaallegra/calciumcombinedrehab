clc
matlabFolder = '/Users/alessandro/Desktop/ELABORAZIONE DATA/';
startFolder = [matlabFolder, '05_Maximum_Intensity_Projection_SeqSelez/_Store_MIP_SIP_all_Days_and_MAT/'];
subjectList = {'GCaMPChR2_30_stroke_Rehab', 'GCaMPChR2_31_stroke_Rehab', 'GCaMPChR2_32_stroke_Rehab',...
     };
 
 
% some group information
GroupName='GlobalSiP-ROBOT-NO-FRICTION';
FieldOfView = 4.4;
NumberOfPixels = 512;


weekList = {'1','2','3','4'};

globalSip = [];
for subject = subjectList
    disp(['Subject: ',subject{1}])
    for week = weekList
        disp(['    Week: ',week{1}])
        if week{1} == '_'
            sipFile = ['Im_AlongDays','_',subject{1},'_Long_SelSeq_SIP.mat'];
        else
            
            sipFile = ['Im_AlongDays_Week_',week{1},'_',subject{1},'_Long_SelSeq_SIP.mat'];
        end
            disp(sipFile)
            sipImg = load([startFolder,subject{1},filesep,sipFile]);
            sipImg = sipImg.Im_med;
            if week{1} == '1'
                %             disp('changing')
                %             sipImg(1:40,1:40)=5;
            end
        
        globalSip = cat(3, globalSip, sipImg);
        
            
    end
end


week1Sip = globalSip(:,:,1:4:end);
week2Sip = globalSip(:,:,2:4:end);
week3Sip = globalSip(:,:,3:4:end);
week4Sip = globalSip(:,:,4:4:end);



fig = figure('name', GroupName);
subplot(2,4,1)
imshow(mean(week1Sip,3)/max(max(mean(week1Sip,3))), 'colormap', jet)

title('Week1')
subplot(2,4,2)

imshow(median(week2Sip,3)/max(max(median(week2Sip,3))), 'colormap', jet)
title('Week2')
subplot(2,4,3)
imshow(median(week3Sip,3)/max(max(median(week3Sip,3))), 'colormap', jet)
title('Week3')
subplot(2,4,4)

imshow(median(week4Sip,3)/max(max(median(week4Sip,3))), 'colormap', jet)
title('Week4')


subplot(2,4,5)
toPlot=week1Sip;
toPlot(toPlot<3/5)=0;
toPlot=median(toPlot,3);
toPlot(toPlot<3/5) = 0;
imshow(toPlot, 'colormap', jet)

[x,y] = find(toPlot~=0);
centroid = [mean(x), mean(y)];
weights = toPlot(toPlot~=0);
weights = weights/sum(weights)*size(weights,1);
centroid_w = [mean(x.*weights), mean(y.*weights)];
hold on;
% scatter(centroid(2),centroid(1),300,'filled','MarkerFaceColor','white','Marker','s')
scatter(centroid_w(2),centroid_w(1),100,'filled','MarkerFaceColor','white','MarkerEdgeColor','black','Marker','o')
hold off;
bdist=sqrt((1-centroid_w(2))^2+(29-centroid_w(1))^2)*FieldOfView/NumberOfPixels;
area = (FieldOfView/NumberOfPixels)^2 * size(x,1);
title({'Week1','dist mm, Area mm^2',[num2str(bdist),',',num2str(area)]})

xlabel(fig.Name)

h=fig;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,[fig.Name,'_averageActivation'],'-dpdf','-r0')
subplot(2,4,6)
toPlot=week2Sip;
toPlot(toPlot<3)=0;
toPlot=median(toPlot,3)/max(max(median(toPlot,3)));
toPlot(toPlot<.60) = 0;
imshow(toPlot, 'colormap', jet)
title('Week2')
[x,y] = find(toPlot~=0);
centroid = [mean(x), mean(y)];
weights = toPlot(toPlot~=0);
weights = weights/sum(weights)*size(weights,1);
centroid_w = [mean(x.*weights), mean(y.*weights)];
hold on;
scatter(centroid(2),centroid(1),300,'filled','MarkerFaceColor','white','Marker','s')
scatter(centroid_w(2),centroid_w(1),150,'filled','MarkerFaceColor','green','Marker','d')
hold off;
bdist=sqrt((1-centroid_w(2))^2+(29-centroid_w(1))^2)*FieldOfView/NumberOfPixels;
area = (FieldOfView/NumberOfPixels)^2 * size(x,1);
title({'Week2','dist mm, Area mm^2',[num2str(bdist),',',num2str(area)]})

subplot(2,4,7)
toPlot=week3Sip;
toPlot(toPlot<3)=0;
toPlot=median(toPlot,3)/max(max(median(toPlot,3)));
toPlot(toPlot<.60) = 0;
imshow(toPlot, 'colormap', jet)
title('Week3')
[x,y] = find(toPlot~=0);
centroid = [mean(x), mean(y)];
weights = toPlot(toPlot~=0);
weights = weights/sum(weights)*size(weights,1);
centroid_w = [mean(x.*weights), mean(y.*weights)];
hold on;
scatter(centroid(2),centroid(1),300,'filled','MarkerFaceColor','white','Marker','s')
scatter(centroid_w(2),centroid_w(1),150,'filled','MarkerFaceColor','green','Marker','d')
hold off;
bdist=sqrt((1-centroid_w(2))^2+(29-centroid_w(1))^2)*FieldOfView/NumberOfPixels;
area = (FieldOfView/NumberOfPixels)^2 * size(x,1);
title({'Week3','dist mm, Area mm^2',[num2str(bdist),',',num2str(area)]})

subplot(2,4,8)
toPlot=week4Sip;
toPlot(toPlot<3)=0;
toPlot=median(toPlot,3)/max(max(median(toPlot,3)));
toPlot(toPlot<.60) = 0;
imshow(toPlot, 'colormap', jet)
[x,y] = find(toPlot~=0);
centroid = [mean(x), mean(y)];
weights = toPlot(toPlot~=0);
weights = weights/sum(weights)*size(weights,1);
centroid_w = [mean(x.*weights), mean(y.*weights)];
hold on;
scatter(centroid(2),centroid(1),300,'filled','MarkerFaceColor','white','Marker','s')
scatter(centroid_w(2),centroid_w(1),150,'filled','MarkerFaceColor','green','Marker','d')
scatter(1,29)
hold off;
bdist=sqrt((1-centroid_w(2))^2+(29-centroid_w(1))^2)*FieldOfView/NumberOfPixels;
area = (FieldOfView/NumberOfPixels)^2 * size(x,1);
title({'Week4','dist mm, Area mm^2',[num2str(bdist),',',num2str(area)]})


savefig(fig, [fig.Name,'_averageActivation.fig'])


%%
figSub=figure('name', GroupName);
weeki=0;
weekSipList = {week1Sip, week2Sip, week3Sip, week4Sip};
for week= weekSipList
    weeki=weeki+1;
    
    for sub=1:length(subjectList)
        

        toPlot=(week{1}(:,:,sub))/max(max(week{1}(:,:,sub)));
        toPlot(toPlot<3/5)=0;
        subplot(4,6,sub+(weeki-1)*6)
        imshow(toPlot, 'colormap', jet)
        hold on
        [x,y] = find(toPlot~=0);
        
        centroid = [mean(x), mean(y)];
        weights = toPlot(toPlot~=0);
        weights = weights/sum(weights)*size(weights,1);
        centroid_w = [mean(x.*weights), mean(y.*weights)];
        area = sum((FieldOfView/NumberOfPixels)^2 * toPlot(sub2ind(size(toPlot),x,y)));
        bdist=sqrt((1-centroid_w(2))^2+(29-centroid_w(1))^2)*FieldOfView/NumberOfPixels;
        
        
        if sub == 1
            ylabel(['Week ',num2str(weeki)])
        end
        if weeki == length(weekSipList)
            xlabel(['Subject ', subjectList{sub}])
        end
        if weeki == 1
            title([subjectList{sub}],'interpreter','None')
        end
        scatter(centroid(2),centroid(1),300,'filled','MarkerFaceColor','white','Marker','s')
        scatter(centroid_w(2),centroid_w(1),150,'filled','MarkerFaceColor','green','Marker','d')
        title({'Dist mm, Area mm^2',[num2str(bdist),',',num2str(area)]})
    end
    
end
% 
% %%
% savefig(fig, 'averageActivation.fig')
% savefig(figSub, 'activationBySubject.fig')