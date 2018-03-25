function res=extractPos(filename,videopath,threshold,Px,PLED,...
                        rect_Marker, rect_LED)

filtwind=6;%15
framenumber=1;


 obj = mmreader(filename);
 totalframes = obj.NumberOfFrames;
 
 if( framenumber <=  totalframes)     
 	video = read(obj,framenumber);
    
        
    grayImm = rgb2gray(video(...
                    rect_Marker(2):rect_Marker(2)+rect_Marker(4),...
                    rect_Marker(1):rect_Marker(1)+rect_Marker(3),...
                    :));
    
   
                
    grayImmLED = rgb2gray(video(...
                    rect_LED(2):rect_LED(2)+rect_LED(4),...
                    rect_LED(1):rect_LED(1)+rect_LED(3),...
                    :));            
    %
    filtImm = medfilt2(grayImm,[filtwind filtwind]);
    filtImmLED = medfilt2(grayImmLED,[filtwind filtwind]);
    
    LEDgrayscale(1) = filtImmLED(round(PLED(2)),round(PLED(1)));
    
    bw = im2bw(filtImm,threshold);
    
    se = strel('disk',4);
    imdil = imdilate(bw,se);
    imdiler = imerode(imdil,se);
    bw = imdiler;
 end
 
%N.b. 1 = bianco; 0 = nero 
% % lenx = 20;
% % leny = 20;


X = zeros(totalframes,2);
LED = zeros(totalframes,1);

X(1,:) = Px + [rect_Marker(1) rect_Marker(2)];
X_buf(1,:) = Px;
LED(1)=0;

%visualizza waitbar
h = waitbar(0,'Please wait...');

for i=2:totalframes
    
    %
    lenx = 20;
    leny = 20;
    %
    
    display(strcat(num2str((i/totalframes)*100),'_% completed'))
    framenumber=i
    % aggiorna waitbar 
    waitbar(framenumber / totalframes)
    
    video = read(obj,framenumber);
    
%     grayImm = rgb2gray(video);
    %
    grayImm = rgb2gray(video(...
                    rect_Marker(2):rect_Marker(2)+rect_Marker(4),...
                    rect_Marker(1):rect_Marker(1)+rect_Marker(3),...
                    :));
                
    grayImmLED = rgb2gray(video(...
                    rect_LED(2):rect_LED(2)+rect_LED(4),...
                    rect_LED(1):rect_LED(1)+rect_LED(3),...
                    :));            
    %
    filtImm = medfilt2(grayImm,[filtwind filtwind]);
    filtImmLED = medfilt2(grayImmLED,[filtwind filtwind]);
    
    LEDgrayscale(i) = filtImmLED(round(PLED(2)),round(PLED(1)));
    
    bw = im2bw(filtImm,threshold);
    
    se = strel('disk',4);
    imdil = imdilate(bw,se);
    imdiler = imerode(imdil,se);
    bw = imdiler;

%     imshow(bw)
    % nuovo centroide della vecchia finestra Mx
    Mx = bw(round(Px(2))-leny:round(Px(2))+leny, round(Px(1))-lenx:round(Px(1))+lenx);
    [bound,label] = bwboundaries(Mx);
%     label
    
    %nel caso in cui il marker nero avesse suoi valori nel bordo della
    %finestra, non si riesce a estrarre il centroide. Allora viene
    %allargata la finestra lungo y per far cadere l'intero marker nero
    %all'interno della finestra
%if (i~=114) && (i
%       while   logical( length (find(label == 0) ))   &&     (~ logical( length (find(label >= 2) )))
      if logical( length (find(label == 0) ))   &&     (~ logical( length (find(label >= 2) )))
        lenx=lenx+5;
        leny=leny+5;
        Mx = bw(round(Px(2))-leny:round(Px(2))+leny, round(Px(1))-lenx:round(Px(1))+lenx);
        [bound,label] = bwboundaries(Mx);
%         label
      end
            
    %
    
    %label è un array 2D di interi non-negativi che rappresenta le regioni
    %continue. Di questa regione si cerca il centroide
    temp1 = regionprops(label,'Centroid');
    
    %%%%nel caso il centroide non venga trovato utilizza la posizione
    %%%%precedente
    if size(temp1,1)>1
        res1 = temp1(2,1).Centroid; %ok
        deltax = res1(1)-lenx;      %ok
        deltay = res1(2)-leny;      %ok
    else
        deltax = 0;
        deltay = 0;        
    end
    %%%
    Px = X_buf(i-1,:)+[deltax deltay];
    X_buf(i,:)=Px;

%     pause(0.2)
    
    X(i,:)=Px + [rect_Marker(1) rect_Marker(2)]
%     else
%         X_buf(i,:)=Px;
%         X(i,:)=Px;
%     end
        
    
end

%chiude waitbar
close(h);

X_temp = zeros(size(X));
X_temp(:,1) = X(:,2);
X_temp(:,2) = X(:,1);

figure(1),
subplot(2,1,1), 
plot(0:1/25:(totalframes-1)/25,X_temp(:,1))
subplot(2,1,2), 
plot(0:1/25:(totalframes-1)/25,X_temp(:,2))

res.X = X_temp;
res.LEDgrayscale= LEDgrayscale;

CurrDir = pwd;
cd(videopath)
save(strcat(videopath,'Pos.mat'),'res')
cd(CurrDir)
