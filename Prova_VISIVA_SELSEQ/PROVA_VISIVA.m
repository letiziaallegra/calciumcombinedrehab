% Script Matlab per selezionare e scartare le sequenze di frames centrate
% intorno a un picco di forza
%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

%%%
Animale_Start  = 'GCampChR2_TOX1';
% Animale_Start  = [];

%%%
CURRFOLDER     = cd;
%WORKING_FOLDER = [CURRFOLDER,'\','Sequenze_Trials_DATA'];
WORKING_FOLDER = ['C:\Users\asus\Desktop\Sequence_Trials_DATA'];
%%%
% Day_Start      = 'DAY';


%%%%%%%%%%%%%%%%%%%%%%%

all_folder_animals = dir(WORKING_FOLDER);

for i_fA = 3:length(all_folder_animals) %for all animals
    
    animal_Name   = all_folder_animals(i_fA,1).name;
    Len_Animal_Name = length(animal_Name);
    
%     if (strfind(animal_Name,'GCaMP') && isempty(Animale_Start)) || strcmp(animal_Name,Animale_Start) %if animal
    if (contains(animal_Name,'GCamp') && isempty(Animale_Start)) || strcmp(animal_Name,Animale_Start) %if animal
        
        all_folder_day = dir([WORKING_FOLDER,'\',animal_Name]);
        
        for i_fD = 3:length(all_folder_day) %for all days
            
            folder_day = all_folder_day(i_fD,1).name;
            
            %if (strfind(folder_day,'GCaMP')) %if day
            if (strfind(folder_day,'GCamp')) %if day 
                %%% day
                dayNum = folder_day(Len_Animal_Name+2:Len_Animal_Name+3);
                
                %%% current folder
                Curr_folder = [WORKING_FOLDER,'\',animal_Name,'\',animal_Name,'_',dayNum,'_SequenceLong_TIF','\SequenceLong'];
                %%% before loading ImageSequence
                filenameDataSeq = [Curr_folder,'\','ImageSequence_',animal_Name,'_',dayNum,'_CenFrame_1_3_ROI_1'];
                
                if exist([filenameDataSeq,'.mat'],'file') == 2
                    %%% loading ImageSequence
                    load(filenameDataSeq)
                    
                    %%%%%%%%% ANALYSIS %%%%%%%%
                    
                    x_all         = zeros(size(ImageSequence.MatrixImageSequence,1),1);
                    y_all         = zeros(size(ImageSequence.MatrixImageSequence,1),1);
                    Ok_No_par_all = zeros(size(ImageSequence.MatrixImageSequence,1),1);
                    
                    for i_imS=1:size(ImageSequence.MatrixImageSequence,1)
                        
                        ImageStack = ImageSequence.MatrixImageSequence{i_imS,1};
                        NumImage   = size(ImageStack,3);
                        
                        figSubTrial = figure('Name',['Trial n. ', num2str(i_imS), ' of ', num2str(size(ImageSequence.MatrixImageSequence,1)), ' (',animal_Name,'_',dayNum,')']);
                        fr          = 1;
                        index_plot  = 1;
                        
                        MAX_plot = max(max(max(ImageStack)));
                        MIN_plot = min(min(min(ImageStack)));
                        
                        for j=1:NumImage
                            
                            if j == fr
                                Image = ImageStack(:,:,j);
                                
%                                 if j==1
%                                     Image_FirstFrame = Image;
%                                 else
%                                     
%                                     Image = Image - Image_FirstFrame;
%                                 end
                                
                                subplot(5,4,index_plot)
                                imagesc(Image)
                                colorbar
                                %caxis([MIN_plot MAX_plot])
                                caxis([-50 50])
                                
                                index_plot = index_plot+1;
                                fr = fr+3;
                            end
                            
                        end
                        [x y] = ginput;
                        
                        close
                        
                        if ~isempty(x)
                            x = round(x(end));
                            y = round(y(end));
                            
                            %if you have selected the image centroid you want to keep it
                            Ok_No_par = 1;
                            
                        else
                            x = 0;
                            y = 0;
                            
                            %if it is empty you want to trow it away
                            Ok_No_par = 0;
                        end
                        
                        x_all(i_imS)           = x;
                        y_all(i_imS)           = y;
                        Ok_No_par_all(i_imS)   = Ok_No_par;
                        
                    end
                    %save DataInfoSequence
                    DataInfoSequence = [x_all,  y_all, Ok_No_par_all];
                    filename_Data = ['DataInfoSequence_', animal_Name, '_', dayNum];
                    save([Curr_folder,'\',filename_Data], 'DataInfoSequence');
                    
                    %save DataInfoSequence
                    figure('Name','Scatter Plot')
                    scatter(y_all,x_all);
                    filename_Image = ['Scatter Plot_',animal_Name, '_', dayNum,'_FIG'];
                    saveas(gca,[Curr_folder,'\',filename_Image],'fig');
                    saveas(gca,[Curr_folder,'\',filename_Image],'jpeg');
                    close
                    
                     display(['Animal: ', animal_Name,' Day: ', dayNum,' done and saved']);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                else
                    %%% file missing
                    display(['Animal: ', animal_Name,' Day: ', dayNum,' is missing']);
                end
                
            end %end if day
            
        end %end for all days
        
    end %end if animals
    
end %end all animals

display('End of the process')

    