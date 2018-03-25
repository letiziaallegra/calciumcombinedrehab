% function res = batch(videofile,datafile)
function res = batch

DirCurrent = cd;

[videofile videopath] = ...
        uigetfile({'*.avi','Video Files (*.avi)'; ...
        '*.*', 'All Files (*.*)'}, 'Select the Video File');

cd(videopath)
[datafile datapath] = ...
        uigetfile({'*.txt','Data Files (*.txt)'; ...
        '*.*', 'All Files (*.*)'}, 'Select the Data File');
    
targetPos = menu('Target Position','8mm','10mm');
targetPos = (targetPos==1)*8 + (targetPos==2)*10;
    
cd(DirCurrent)

% [bw BWthres] = setBWthreshold(videofile);
MaxGrayLevel= menu('MaxGrayLevel to use','100: Morning or Sunny Day','75: Afternoon or ');
MaxGrayLevel = (MaxGrayLevel==1)*100 + (MaxGrayLevel==2)*75;
% MaxGrayLevel = 100;
[bw BWthres] = setBWthreshold_v2([videopath videofile],MaxGrayLevel);

display('Threshold set')

% G = viewSingleFrame(videofile,1,'','gray');
G = viewSingleFrame([videopath videofile],1,'','gray');

[PLED, P, rect_Marker, rect_LED] = Snap(bw,G);
display('Marker e LED images extracted')

 %Tolto il 18 maggio
% P = findPositionMarker(bw);
% display('Marker position found')

% resPos = extractPos_v2([videopath videofile],videopath,BWthres,P,PLED,...
%                     rect_Marker, rect_LED);
                
resPos = extractPos_v2([videopath videofile],videopath,BWthres,P,PLED,...
                    rect_Marker, rect_LED);

display('Position extracted')
% resSync = synchronize('Pos.mat',datafile);
[resSync newfilename]= synchronize('Pos.mat',[datapath datafile], datapath,targetPos);
display('Signals Synced')
cd(datapath)

%plot della posizione estratta
txt = load(newfilename,'ASCII');
figure
plot(txt(:,1),txt(:,7)/20,'b')
hold on
plot(txt(:,1),txt(:,8),'g')
hold on
clear txt datamouse newfilename

