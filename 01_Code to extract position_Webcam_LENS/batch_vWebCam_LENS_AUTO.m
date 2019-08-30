%%% batch_VIDEOCAMERA_AUTO
clear
clc

FunDir  = pwd;

%PATH
%DirCurr = 'M:\LENS\Animals Data STIM\170511\RCAMP1A 11.04';
%DirCurr = '/Users/alessandro/DATA/Imaging/matlab_sample_data/GCaMP24/';
DirCurr = '/Volumes/ALE6TB_DESK/DATA/Imaging/Emilia/toxin/MATLAB/Tox1_toxin';
DirCurr = '/Users/alessandro/Desktop/180424_RehabOptogen/MATLAB/or27_robot';
DirCurr = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB/';
DirCurr = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_GCaMP27-29';
DirCurr = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_OR';
DirCurr = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_CTRL_190122';
DirCurr = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_CTRL_190117';

a_dir = dir(DirCurr);
for a_sub_dir=3:length(a_dir)
    disp(a_dir(a_sub_dir).name)
    ListSubFolder = dir([DirCurr, filesep, a_dir(a_sub_dir).name]);
    
    for i=3:length(ListSubFolder)
        
        close all
        
        indexDay = i-2;
        
        NameDay = ListSubFolder(i,1).name;
        disp(NameDay)
        PathFile = [DirCurr,filesep,a_dir(a_sub_dir).name,filesep,NameDay];
        
        ListFile = dir(PathFile);
        
        if length(ListFile)>2 %not empty
            
            videofile = [];
            datafile  = [];
            
            for j=3:length(ListFile)
                
                Namefile = ListFile(j,1).name;
                if length(Namefile)>6
                    if strcmp(Namefile(end-2:end),'avi') || strcmp(Namefile(end-2:end),'mpg')
                        
                        videofile = Namefile;
                        
                    elseif strcmp(Namefile(end-2:end),'txt') || strcmp(Namefile(end-2:end),'dat')
                        
                        datafile = Namefile;
                    end
                    
                    if j == length(ListFile)
                        if ~isempty(videofile) && ~isempty(datafile)
                            
                            batch_vWebCam_LENS_fun(PathFile,videofile, datafile)
                            cd(FunDir)
                            
                        else
                            display(['No video or No text in ',PathFile])
                        end
                    end
                    
                end
            end
            
            
        else
            display([PathFile, 'is empty'])
        end
    end
end



