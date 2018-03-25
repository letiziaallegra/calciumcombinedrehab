%%% batch_VIDEOCAMERA_AUTO
FunDir  = pwd;
%PATH
DirCurr = 'M:\LENS\Animals Data STIM\170511\RCAMP1A 11.04';
DirCurr = '/Users/alessandro/DATA/Imaging/matlab_sample_data/GCaMP24/';

ListSubFolder = dir(DirCurr);

for i=3:length(ListSubFolder)
    
    close all 
    
    indexDay = i-2;
    
    NameDay = ListSubFolder(i,1).name;
    
    PathFile = [DirCurr,filesep,NameDay];
    
    ListFile = dir(PathFile);
    
    if length(ListFile)>2 %not empty
        
        videofile = [];
        datafile  = [];
        
        for j=3:length(ListFile)
            
            Namefile = ListFile(j,1).name;
            
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
        
                
    else
        display([PathFile, 'is empty'])
    end
end
    



