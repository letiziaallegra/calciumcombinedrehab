function batch_vWebCam_LENS_fun(PathFile,videofile, datafile)

% % DirCurrent = cd;
% % 
% % [videofile videopath] = ...
% %         uigetfile({'*.avi','Video Files (*.avi)'; ...
% %         '*.*', 'All Files (*.*)'}, 'Select the Video File');
% % 
% % cd(videopath)
% % [datafile datapath] = ...
% %         uigetfile({'*.txt','Data Files (*.txt)'; ...
% %         '*.*', 'All Files (*.*)'}, 'Select the Data File');
% %     
% %     
% % cd(DirCurrent)

datapath  = PathFile;
videopath = PathFile;

obj = VideoReader([videopath,filesep,videofile]); 
frame = read(obj,1);
grayImm = rgb2gray(frame); 
I = medfilt2(grayImm,[6 6]);

%marker area and LED detection
[P_Marker, rect_Marker, rect_LED] = Snap_vWebCam_LENS_auto(I);
display('Marker e LED images extracted')

%position extraction
resPos = extractPos_vWebCam_LENS([videopath,filesep,videofile],videopath,P_Marker,rect_Marker, rect_LED);
display('Position extracted')

% resSync = synchronize('Pos.mat',datafile);
targetPos = 8 ;
[resSync newfilename]= synchronize_vWebCam_LENS_v3('Pos.mat',[datapath,filesep,datafile], datapath,targetPos);
display('Signals Synced')
cd(datapath)

%plot della posizione estratta
txt = load(newfilename,'ASCII');
figure('Name',newfilename)
plot(txt(:,1),txt(:,7),'b')
hold on
plot(txt(:,1),txt(:,9),'g')
hold on
saveas(gca,strcat(videopath,'\SynchroPositionFig'),'jpeg')
clear txt datamouse newfilename



end
