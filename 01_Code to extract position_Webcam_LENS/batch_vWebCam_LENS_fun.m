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

% % added fast skip if video already analyzed with python algs
python_folder_results_file = strcat(videopath,filesep,'analysis_behavior',filesep,'position_and_force.csv');
if exist(python_folder_results_file,'file')
    disp(['found python results'])
    % loading results
    table = readtable(python_folder_results_file);
    data = table2array(table(:,1:9));
    % the following code has been copied from synchronize_vWebCam_LENS_v3
    % finding time 0 start of the fluorescence
    find(data(:,1) == 0)
    Pos_filt = data(:,9);
    Vel = derivative(Pos_filt,0.01);
    Vel_filt = sgolayfilt(Vel,3,11);
    Acc = derivative(Vel_filt,0.01);
    Acc_filt = sgolayfilt(Acc,3,11);
    filenameFORCE_RAW = [datapath,filesep,datafile];
    newfilename = strcat(filenameFORCE_RAW(1:end-4),'_sync.txt');
    
    t_res = 0:0.01:(length(Pos_filt)-1)*0.01;
    res.t_res = t_res';
    res.Pos=Pos_filt;
    res.Vel=Vel_filt;
    res.Acc=Acc_filt;
    
    data(:,9:11) = [res.Pos res.Vel res.Acc];
    data = data(1:end-1,:);
    
    if(fopen(newfilename)==-1)
        save(newfilename, 'data', '-ascii');
    else
        fclose('all')
        delete(newfilename);
        save(newfilename, 'data', '-ascii');
    end
    
    cd(datapath)

    %plot della posizione estratta
    txt = load(newfilename,'ASCII');
    figure('Name',newfilename)
    plot(txt(:,1),txt(:,7),'b')
    hold on
    plot(txt(:,1),txt(:,9),'g')
    hold on
    saveas(gca,strcat(videopath,filesep,'SynchroPositionFig'),'jpeg')
    clear txt datamouse newfilename
    
    return
end


sync_file = strcat(videopath,filesep,'SynchroPositionFig.jpg');

if exist(sync_file) == 2
    disp(['Synch File found, skipping:',datapath])
    return
end

obj = VideoReader([videopath,filesep,videofile]);
frame = read(obj,1);
grayImm = rgb2gray(frame);
I = medfilt2(grayImm,[6 6]);

%marker area and LED detection
[P_Marker, rect_Marker, rect_LED] = Snap_vWebCam_LENS_auto(I);
display('Marker e LED images extracted')

%position extraction
pos_file = strcat(videopath,filesep,'Pos.mat');
rect_LED = [368 245 270 203];
rect_Marker = [7 152 352 200];
if exist(pos_file)==2
    data = load(pos_file);
    resPos = data.res;
    display(['Position laoded from file. Delete',pos_file,' if you want to recompute the position'])
else
    resPos = extractPos_vWebCam_LENS([videopath,filesep,videofile],videopath,P_Marker,rect_Marker, rect_LED);
    display('Position extracted')
end

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
saveas(gca,strcat(videopath,filesep,'SynchroPositionFig'),'jpeg')
clear txt datamouse newfilename
