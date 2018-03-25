%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Costruisci sequenza di Im, a 8 bits e Flippata Vertically %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%le figure sono salvate come matrici .mat, 1 per Immagine in una cartella 

clear
close all
clc

%MainDir       = 'I:\LENS\Script_FlipStacks_Auto\Animals_Data_Make_Mat_Seq_to_do';
MainDir       = 'C:\LENS\Data Leti';
MainDirFolder = dir(MainDir);

for mdf=3:length(MainDirFolder) %for MainDirFolder
    
    AnimalDir       = MainDirFolder(mdf,1).name;
    AnimalDirFolder = dir([MainDir,'\',AnimalDir]);
    
    for adf=3:length(AnimalDirFolder)  %for AnimalDirFolder
        
        DayDir       = AnimalDirFolder(adf,1).name;
        DayDirFolder = dir([MainDir,'\',AnimalDir,'\',DayDir]);
        
        for adf=3:length(DayDirFolder) %for DayDirFolder
            
            FileDayDir   = DayDirFolder(adf,1).name;
            
            if strcmp(FileDayDir(1:5),'trial') %if FileDayDir
                
                %check the image dir
                if isdir([MainDir,'\',AnimalDir,'\',DayDir,'\',FileDayDir]); %if check the image dir
                    
                    ImageDir = dir([MainDir,'\',AnimalDir,'\',DayDir,'\',FileDayDir]);                    
                    
                    %check if the dir is full
                    if length(ImageDir)>3 %if dir is full
                        
                        %to have a sorted list of matlab file, I need to
                        %scroll the images following the right order 
                        %1,2,3 not 1, 1000, 1001 and so on
                        MaxNum_Str = length(num2str(length(ImageDir)-2));
                        
                        for imd=3:length(ImageDir) %for ImageDir
                                                        
                            ImageFileName     = ImageDir(imd,1).name;
                            UndSc_pos         = strfind(ImageFileName ,'_');
                            UndSc_pos         = UndSc_pos(end);
                            Dot_pos           = strfind(ImageFileName ,'.');
                            Dot_pos           =  Dot_pos(end);
                            Num_str           = ImageFileName(UndSc_pos+1:Dot_pos-1);
                            ImageFileNameRoot = ImageFileName(1:UndSc_pos);                    
                            
                            Im16 = importdata([MainDir,'\',AnimalDir,'\',DayDir,'\',FileDayDir,'\',ImageFileName]);
                            
                            if ~iscell(Im16)
                                
                                %16 bits -> 8 bits
                                M16 = 2^16-1;
                                M8 =  2^8-1;
                                Im8 = uint8(Im16 * (M8/M16));
                                
                                %flip vertically
                                Im8_fv = flipdim(Im8 ,1);
                                
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
                                Folder_MAT_Name = [MainDir,'\',AnimalDir,'\',DayDir,'\','MAT_trial'];
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
                                
                                save([Folder_MAT_Name,'\',ImageFileNameMAT],'Im8_fv')
                            else
                                
                                display([AnimalDir,'\',DayDir,'\',FileDayDir,'\',ImageFileName, ' is a corrupted image-file'])
                            end
                            
                            
                            
                        end %end ImageDir
                        
                        break
     
                    else
                        error([MainDir,'\',AnimalDir,'\',DayDir,'\',FileDayDir,' is empty'])
                        
                    end %end if dir is full
                    
                else 
                    display([MainDir,'\',AnimalDir,'\',DayDir,'\',FileDayDir,' is not a directory'])
                    
                end %end if check the image dir
                                   
            end %end if FileDayDir
            
        end %end for DayDirFolder
        display(['END DAY ', DayDir,' of ANIMAL ', AnimalDir])
        
    end %end for AnimalDirFolder
    display(['END ANIMAL ',  AnimalDir])
    
end  %end for MainDirFolder

display('END PROCESS')