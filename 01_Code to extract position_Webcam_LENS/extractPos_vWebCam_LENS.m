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
fh = figure();
video = rgb2gray(readFrame(obj));

video(rect_Marker(2):rect_Marker(2)+rect_Marker(4),rect_Marker(1))=256;
video(rect_Marker(2),rect_Marker(1):rect_Marker(1)+rect_Marker(3))=256;
video(rect_Marker(2):rect_Marker(2)+rect_Marker(4),rect_Marker(1)+rect_Marker(3))=256;
video(rect_Marker(2)+rect_Marker(4),rect_Marker(1):rect_Marker(1)+rect_Marker(3))=256;

video(rect_LED(2):rect_LED(2)+rect_LED(4),rect_LED(1))=256;
video(rect_LED(2),rect_LED(1):rect_LED(1)+rect_LED(3))=256;
video(rect_LED(2):rect_LED(2)+rect_LED(4),rect_LED(1)+rect_LED(3))=256;
video(rect_LED(2)+rect_LED(4),rect_LED(1):rect_LED(1)+rect_LED(3))=256;
imshow(video)




for i=2:totalframes
    obj.set('CurrentTime',0)
    
    framenumber=i;
    % aggiorna waitbar 
    waitbar(framenumber / totalframes)
    
    %load frame
    obj.CurrentTime=framenumber/obj.FrameRate;
    video = rgb2gray(readFrame(obj));
    %video = read(obj,framenumber);

    %rgb->gray
    grayImm = video(...
                    rect_Marker(2):rect_Marker(2)+rect_Marker(4),...
                    rect_Marker(1):rect_Marker(1)+rect_Marker(3),...
                    :);
                
    grayImmLED = video(...
                    rect_LED(2):rect_LED(2)+rect_LED(4),...
                    rect_LED(1):rect_LED(1)+rect_LED(3),...
                    :);    
    if mod(i,10)==0         
        display(strcat(num2str((i/totalframes)*100),['_% completed ',videopath]))
        video(rect_Marker(2):rect_Marker(2)+rect_Marker(4),rect_Marker(1))=256;
        video(rect_Marker(2),rect_Marker(1):rect_Marker(1)+rect_Marker(3))=256;
        video(rect_Marker(2):rect_Marker(2)+rect_Marker(4),rect_Marker(1)+rect_Marker(3))=256;
        video(rect_Marker(2)+rect_Marker(4),rect_Marker(1):rect_Marker(1)+rect_Marker(3))=256;

        video(rect_LED(2):rect_LED(2)+rect_LED(4),rect_LED(1))=256;
        video(rect_LED(2),rect_LED(1):rect_LED(1)+rect_LED(3))=256;
        video(rect_LED(2):rect_LED(2)+rect_LED(4),rect_LED(1)+rect_LED(3))=256;
        video(rect_LED(2)+rect_LED(4),rect_LED(1):rect_LED(1)+rect_LED(3))=256;
        imshow(video)
    end
    %filtering
    % filtImm = medfilt2(grayImm,[filtwind filtwind]);
    filtImmLED = medfilt2(grayImmLED,[filtwind filtwind]);
     
    %bw Marker
    %bw = im2bw(grayImm,245/255);
    bw = im2bw(grayImm,20/255);
    
    %erode Marker
    %se = strel('disk',10); 
    se = strel('disk',18);
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
try
    close(h)
catch
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

%%fixes the different frame rate
fs = 25; %assumed to be 25
if obj.FrameRate ~= fs
    ra = rats(fs/obj.FrameRate);
    pq=double(split(ra,'/'));
    if pq ~= 1
        LED = resample(LED,pq(1),pq(2),0);
        X_temp = resample(X_temp,pq(1),pq(2),0);
    end
end

res.X = X_temp;
res.LEDgrayscale= LED;




%some data cleaning Alessandro Scaglione 180329
%clearly the LED cannot be on for small timewindows when the position is
%not at maximum or minimum. Removing those data points if any.
intervals = reshape(find(diff(LED)~=0),2,[])';
durations = diff(intervals,1,2); % I may not use this at the end
%for now takes care only of very shor intervals
for sd = 1:size(durations,1)
    temp_int = intervals(sd,:);
    d_length(sd) = diff(temp_int);
    pos_level(sd) = mean(res.X(temp_int(1):temp_int(2),1));
    if d_length(sd) <= 3
        res.LEDgrayscale(temp_int(1):temp_int(2)+1) = 0;
    end
end


    


CurrDir = pwd;
cd(videopath)
if strcmp(videopath(end),filesep)
    save(strcat(videopath,'Pos.mat'),'res')
else
    save(strcat(videopath,filesep,'Pos.mat'),'res')
end
cd(CurrDir)


