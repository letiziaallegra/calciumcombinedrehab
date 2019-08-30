%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Costruisci sequenza di Im, a 8 bits e Flippata Vertically %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%le figure sono salvate come matrici .mat, 1 per Immagine in una cartella

clear
close all
clc

%MainDir       = 'I:\LENS\Script_FlipStacks_Auto\Animals_Data_Make_Mat_Seq_to_do';
%MainDir       = 'C:\LENS\Data Leti';
%MainDir = '/Users/alessandro/DATA/Imaging/matlab_sample_data';
MainDir = '/Volumes/ALE6TB_DESK/DATA/Imaging/Emilia/Rehab/MATLAB_DATA_FOLDERS';
MainDir = '/Volumes/ALE6TB_DESK/DATA/Imaging/Emilia/toxin/MATLAB';
MainDir = '/Users/alessandro/Desktop/180424_RehabOptogen/MATLAB';
MainDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB';
MainDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_GCaMP27-29';
MainDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_OR';
MainDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_CTRL_190122'; FLIP = true;
MainDir = '/Volumes/ALE6TB_DESK/Projects/Imaging/DATA/OptoStimRehabFull/MATLAB_CTRL_190117'; FLIP = true;
MainDirFolder = dir(MainDir);



for mdf=3:length(MainDirFolder) %for MainDirFolder
    
    AnimalDir       = MainDirFolder(mdf,1).name;
    AnimalDirFolder = dir([MainDir,filesep,AnimalDir]);
    
    for adf=3:length(AnimalDirFolder)  %for AnimalDirFolder
        
        DayDir       = AnimalDirFolder(adf,1).name;
        DayDirFolder = dir([MainDir,filesep,AnimalDir,filesep,DayDir]);
        
        for adf=3:length(DayDirFolder) %for DayDirFolder
            
            FileDayDir   = DayDirFolder(adf,1).name;
            if length(FileDayDir)>=3
                %second parts takes out the .analyses folders
                if (strcmp(FileDayDir(1:3),'tri') || strcmp(FileDayDir(1:2),'tx')) && ((length(FileDayDir)<14) || contains(FileDayDir, 'deep')) %if FileDayDir
                    disp(['Importing ',FileDayDir])
                    %check the image dir
                    if isdir([MainDir,filesep,AnimalDir,filesep,DayDir,filesep,FileDayDir]); %if check the image dir
                        
                        ImageDir = dir([MainDir,filesep,AnimalDir,filesep,DayDir,filesep,FileDayDir]);
                        
                        %check if the dir is full
                        if length(ImageDir)>3 %if dir is full
                            
                            %to have a sorted list of matlab file, I need to
                            %scroll the images following the right order
                            %1,2,3 not 1, 1000, 1001 and so on
                            MaxNum_Str = length(num2str(length(ImageDir)-2));
                            
                            
                            Folder_MAT_Name = [MainDir,filesep,AnimalDir,filesep,DayDir,filesep,'MAT_trial'];
                            if ~isdir(Folder_MAT_Name)
                                h = waitbar(0, ['IMPORTING DAY ',DayDir,' of ANIMAL ', AnimalDir]);
                                
                                for imd=3:length(ImageDir) %for ImageDir
                                    
                                    imageName = [ImageDir(imd,1).folder, filesep, ImageDir(imd,1).name];
                                    if ~isdir(imageName)
                                        ImageFileName     = ImageDir(imd,1).name;
                                        UndSc_pos         = strfind(ImageFileName ,'_');
                                        UndSc_pos         = UndSc_pos(end);
                                        Dot_pos           = strfind(ImageFileName ,'.');
                                        Dot_pos           =  Dot_pos(end);
                                        Num_str           = ImageFileName(UndSc_pos+1:Dot_pos-1);
                                        ImageFileNameRoot = ImageFileName(1:UndSc_pos);
                                        %check for file extension to avoid errors
                                        [~,~,ext] = fileparts(ImageDir(imd,1).name);
                                        
                                        if strcmpi(ext,'.tif')
                                            Im16 = importdata([MainDir,filesep,AnimalDir,filesep,DayDir,filesep,FileDayDir,filesep,ImageFileName]);
                                            
                                            if ~iscell(Im16)
                                                
                                                %16 bits -> 8 bits
                                                M16 = 2^16-1;
                                                M8 =  2^8-1;
                                                Im8 = uint8(Im16 * (M8/M16));
                                                
                                                %flip vertically
                                                if FLIP
                                                    Im8_fv = flipdim(Im8 ,1);
                                                else
                                                    Im8_fv = Im8;
                                                end
                                                
                                                %                             % plot
                                                %                             FIG = figure;
                                                %                             subplot(131)
                                                %                             imshow(Im16);
                                                %                             subplot(132)
                                                %                             imshow(Im8);
                                                %                             subplot(133)
                                                %                             imshow(Im8_fv);
                                                %                             iptsetpref('ImshowBorder','tight');
                                                %                             %
                                                
                                                %%% SAVE MAT MATRIX %%%
                                                Folder_MAT_Name = [MainDir,filesep,AnimalDir,filesep,DayDir,filesep,'MAT_trial'];
                                                if ~isdir(Folder_MAT_Name)
                                                    
                                                    %make Folder where saving
                                                    mkdir(Folder_MAT_Name)
                                                end
                                                
                                                %%%%%%  save MATRIX
                                                %%% add zeros to the filename to have (0001) and no (1) (sorted list of image-Matrices)
                                                AddZeroStr = cell(1,[MaxNum_Str - length(Num_str)]);
                                                if ~isempty(AddZeroStr)
                                                    ix=cellfun('isempty',AddZeroStr);
                                                    AddZeroStr(ix)={'0'};
                                                    AddZeroStr = char(AddZeroStr)';
                                                    ImageFileNameMAT = [ImageFileNameRoot,AddZeroStr,Num_str];
                                                else
                                                    ImageFileNameMAT = [ImageFileNameRoot,Num_str];
                                                end
                                                
                                                save([Folder_MAT_Name,filesep,ImageFileNameMAT],'Im8_fv')
                                            else
                                                
                                                display([AnimalDir,filesep,DayDir,filesep,FileDayDir,filesep,ImageFileName, ' is a corrupted image-file'])
                                            end
                                        end
                                        
                                        
                                        waitbar(imd/length(ImageDir))
                                    end
                                end %end ImageDir
                                close(h)
                            end
                            break
                            
                        else
                            error([MainDir,filesep,AnimalDir,filesep,DayDir,filesep,FileDayDir,' is empty'])
                            
                        end %end if dir is full
                        
                    else
                        display([MainDir,filesep,AnimalDir,filesep,DayDir,filesep,FileDayDir,' is not a directory'])
                        
                    end %end if check the image dir
                    
                end %end if FileDayDir
            end
        end %end for DayDirFolder
        display(['END DAY ', DayDir,' of ANIMAL ', AnimalDir])
        
    end %end for AnimalDirFolder
    display(['END ANIMAL ',  AnimalDir])
    
end  %end for MainDirFolder

display('END PROCESS')