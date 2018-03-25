function [bw BWthres] = setBWthreshold_vWebCam(filename,MaxGrayLevel,rect_Marker)
% thres= 0.28
%close all
clc
filtwind=6;%15

obj = mmreader(filename);      
video = read(obj,1);

% grayImm = rgb2gray(frame);  

grayImm = rgb2gray(video(...
                    rect_Marker(2):rect_Marker(2)+rect_Marker(4),...
                    rect_Marker(1):rect_Marker(1)+rect_Marker(3),...
                    :));

filtImm = medfilt2(grayImm,[filtwind filtwind]);
res =filtImm;
imshow(filtImm)
% costruisco l'istogramma
[counts x]=imhist(filtImm);
% ricavo la funzione di distribuzione dei livelli di grigio (continua)
counts_filt = sgolayfilt(counts,3,11);
  
%trovo il massimo della funzione di distribuzione
% [mx imx] = max(counts_filt);
[mx imx] = max(counts_filt(1:200,:)); %solo nella prima parte sotto il 200th level
%calcolo le derivate della funzione di distribuzione
counts_der = derivative(counts_filt,1);
n=[];

%trovo i minimi locali
for i=2:256
    if(counts_der(i-1)<0 && counts_der(i)>0)
        n=[n; i];
    end
end

%dei minimi locali prendo solo quelli alla destra del picco max...
nn= n(find(n>imx));
%... non oltre il livello numero 100
nn_good = nn(find(nn<MaxGrayLevel));
%e di questi scelgo il picco assoluto
[m_ng im_ng] = min(counts_filt(nn_good));
threshold = nn_good(im_ng)/255;
bw = im2bw(filtImm,threshold);
se = strel('disk',4);
imdil = imdilate(bw,se);
imdiler = imerode(imdil,se);
bw = imdiler;
imshow(bw);

im_ngX = [];
while isempty(im_ngX)
    
    ThOK = menu('is Thresholding ok?','Yes','Not');
    
    if ThOK == 2 %NOT   
       plot(counts_filt);
       [im_ngX, im_ngY] = ginput;
       threshold = round(im_ngX)/255;
       bw = im2bw(filtImm,threshold);
       se = strel('disk',4);
       imdil = imdilate(bw,se);
       imdiler = imerode(imdil,se);
       bw = imdiler;
       close
       imshow(bw);
       figure
       imshow(grayImm );
       im_ngX = [];
        

    else %YES
       im_ngX =1;
    end
    
end
BWthres = threshold;
close all

return 
