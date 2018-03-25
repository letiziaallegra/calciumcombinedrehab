function res=extractPos_vWebCam_LENS(filename,videopath, Px, rect_Marker, rect_LED)

filtwind = 6;

obj = VideoReader(filename);
totalframes = round(obj.FrameRate*obj.Duration)-1;

X     = zeros(totalframes,2);
X_buf = zeros(totalframes,2);
LED   = zeros(totalframes,1);

%Absolute position
X(1,:) = Px + [rect_Marker(1) rect_Marker(2)];
%Relative Position
X_buf(1,:) = Px;
%LED swithing on
LED(1)=0;


% %%%CANC
% X(1:6585,1) = Px(1) + rect_Marker(1);
% X(1:6585,2) = Px(2) + rect_Marker(2);
% %Relative Position
% X_buf(1:6585,1) = Px(1) ;
% X_buf(1:6585,2) = Px(2) ;
% %LED swithing on
% LED(1:6585)=0;
% %%%%

%visualizza waitbar
h = waitbar(0,'Please wait...');

for i=2:totalframes
        
    display(strcat(num2str((i/totalframes)*100),'_% completed'))
    framenumber=i;
    % aggiorna waitbar 
    waitbar(framenumber / totalframes)
    
    %load frame
    obj.CurrentTime=framenumber/obj.FrameRate;
    video = readFrame(obj);
    %video = read(obj,framenumber);

    %rgb->gray
    grayImm = rgb2gray(video(...
                    rect_Marker(2):rect_Marker(2)+rect_Marker(4),...
                    rect_Marker(1):rect_Marker(1)+rect_Marker(3),...
                    :));
                
    grayImmLED = rgb2gray(video(...
                    rect_LED(2):rect_LED(2)+rect_LED(4),...
                    rect_LED(1):rect_LED(1)+rect_LED(3),...
                    :));            
    %filtering
    filtImm = medfilt2(grayImm,[filtwind filtwind]);
    filtImmLED = medfilt2(grayImmLED,[filtwind filtwind]);
     
    %bw Marker
    %bw = im2bw(grayImm,245/255);
    bw = im2bw(grayImm,50/255);
    
    %erode Marker
    %se = strel('disk',10); 
    se = strel('disk',20);
    imer = imerode(bw,se);
    
    
     
    %Marker
    [bound,label] = bwboundaries(imer);
    
    if ~isempty(bound)
        Cen = regionprops(label,'Centroid');
        Pos = round([Cen(1,1).Centroid(1) Cen(1,1).Centroid(2)]);
               
        X_buf(i,:)= Pos;
        X(i,:)= Pos + [rect_Marker(1) rect_Marker(2)];
        
    else
        display('marker not found')
        
        X_buf(i,:)= X_buf(i-1,:);
        X(i,:)    = X(i-1,:);
    end
    
    %LED
    [bound_L,label_L] = bwboundaries(filtImmLED);
    
     %if it finds the region
     if ~isempty(bound_L)
             LED(i) = 1;
     else
         LED(i) = 0;
     end     
       
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
res.LEDgrayscale= LED;

CurrDir = pwd;
cd(videopath)
if strcmp(videopath(end),'\')
    save(strcat(videopath,'Pos.mat'),'res')
else
    save(strcat(videopath,'\Pos.mat'),'res')
end
cd(CurrDir)


