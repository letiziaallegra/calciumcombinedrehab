function res=extractPos_vWebCam(filename,videopath,threshold,Px,PLED,...
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
    
    bw = ~im2bw(filtImm,threshold);
    
    se = strel('disk',4);
    imdil = imdilate(bw,se);
    se = strel('disk',2);
    imdiler = imerode(imdil,se);
    bw = ~imdiler;
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
    lenx = 40;
    leny = 30;
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
    
    bw = ~im2bw(filtImm,threshold);
    
    se = strel('disk',4);
    imdil = imdilate(bw,se);
    se = strel('disk',2);
    imdiler = imerode(imdil,se);
    bw = ~imdiler;
%     imshow(bw)
   

    % nuovo centroide della vecchia finestra Mx
    if ((round(Px(2))-leny) >0 )                        &&...
       ((round(Px(2))+leny) < size(bw,1))               &&...
       ((round(Px(1))-lenx) >0 )                        &&...
       ((round(Px(1))+lenx) < size(bw,2))
        
        Mx = bw(round(Px(2))-leny:round(Px(2))+leny, round(Px(1))-lenx:round(Px(1))+lenx);
        [bound,label] = bwboundaries(Mx);
           

        %se non trovo nulla -> updating della soglia (che raddoppia)
        if (logical( length (find(label == 0) ))   &&     (~ logical( length (find(label >= 2) )))) || all(all((label==1)))
           
            upd=2;
            bw = ~im2bw(filtImm,threshold*upd);
            se = strel('disk',4);
            imdil = imdilate(bw,se);
            se = strel('disk',2);
            imdiler = imerode(imdil,se);
            bw = ~imdiler;
            
            Mx = bw(round(Px(2))-leny:round(Px(2))+leny, round(Px(1))-lenx:round(Px(1))+lenx);
            [bound,label] = bwboundaries(Mx);
            
            %se ancora non trovo nulla -> allargo finestra di ricerca
            %nel caso in cui il marker nero avesse suoi valori nel bordo della
            %finestra, non si riesce a estrarre il centroide. Allora viene
            %allargata la finestra lungo y  per far cadere l'intero marker nero
            %all'interno della finestra            
            if (logical( length (find(label == 0) ))   &&     (~ logical( length (find(label >= 2) )))) || all(all((label==1)))
                lenx=lenx+20;
                leny=leny;
                if ((round(Px(2))-leny) >0 )                        &&...
                        ((round(Px(2))+leny) < size(bw,1))          &&...
                        ((round(Px(1))-lenx) >0 )                   &&...
                        ((round(Px(1))+lenx) < size(bw,2))
                    
                    Mx = bw(round(Px(2))-leny:round(Px(2))+leny, round(Px(1))-lenx:round(Px(1))+lenx);
                    [bound,label] = bwboundaries(Mx);
                end
            end
        end



        %label è un array 2D di interi non-negativi che rappresenta le regioni
        %continue. Di questa regione si cerca il centroide
        temp1 = regionprops(label,'Centroid');
    
    else
        temp1 = [];
    end
        


    %%%%nel caso il centroide non venga trovato chiedi all'utente
    if size(temp1,1)>1 
        
        res1 = temp1(2,1).Centroid;         %ok
        deltax = res1(1)-lenx;              %ok
        deltay = res1(2)-leny;              %ok
        
        Px = X_buf(i-1,:)+[deltax deltay];  %ok
        
    else
        %chiede all'utente
        
        viewSingleFrame(filename,framenumber,'plot','gray',threshold);
        %draw rectangular area
        rectangle('Position',[rect_Marker(1),rect_Marker(2), rect_Marker(3),rect_Marker(4)]);
              
        PxUser = [];
        
        while isempty(PxUser)
            
            leny = leny + 5 ;
            lenx = lenx + 5 ;
            rectangle('Position',[...
                                  rect_Marker(1) + round(Px(1)) - lenx,...
                                  rect_Marker(2) + round(Px(2)) - leny,...
                                  2* lenx,...
                                  2* leny...
                                  ]);
                                  
            [PxUser, PyUser] = ginput;
        end
        
        Px = [PxUser(end) - rect_Marker(1),...
              PyUser(end) - rect_Marker(2)];
         
        h = gcf;
        close(h)        
    end
    
    
    X_buf(i,:)=Px;
    X(i,:)= Px + [rect_Marker(1) rect_Marker(2)];

    
end

X_temp = zeros(size(X));
%inverto x con y
X_temp(:,1) = X(:,1);  %in questo caso l'asse x di pulling  è visto dalla telecamera come un asse y 
X_temp(:,2) = X(:,2);  %in questo caso l'asse y di pulling  è visto dalla telecamera come un asse x 

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
