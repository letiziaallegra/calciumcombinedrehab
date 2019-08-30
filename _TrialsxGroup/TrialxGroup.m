clear
clc

AnList = {'GCaMPChR2_7', ...
            'GCaMPChR2_17', ...
            'GCaMPChR2_23', ...
            'GCaMPChR2_24', ...
            ...
            'GCaMPChR2_8',  ...
            'GCaMPChR2_9', ...
            'GCaMPChR2_19', ...
            'GCaMPChR2_22', ...
            'GCaMPChR2_25', ...
            'GCaMPChR2_26', ...
            ...
            'GCaMPChR2_11', ...
            'GCaMPChR2_12', ...
            'GCaMPChR2_14', ...
            'GCaMPChR2_15', ...
            'GCaMP16',       ...
            'GCaMP18'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% %%% %%% %%% %%% %%% %%% %%% %%%
%Folder Data
User = getenv('username');
UsbPortHD  = 'M';
% dirAllDataGCamp = [UsbPortHD,':\LENS\_data_MAT_GCamp\'];
dirAllDataGCamp = ['C:\Users\asus\Google Drive\Piattaforma Stefano\ELABORAZIONE DATA\_data_MAT_GCamp_Store_Analysis_Par'];
dirAllDataGCamp = ['/Users/alessandro/Desktop/ELABORAZIONE DATA/_data_MAT_GCamp_Store_Analysis_Par'];
%list animals
list_dirDataGCamp = dir(dirAllDataGCamp);

%num animal
n_animal_count = 0;

out_cell = {'animal'};
row=2;
for dLGC=3:length(list_dirDataGCamp) %for animals
    
    n_animal_count = n_animal_count+1;
    Treat_old = [];
    
    NameAnimal =  list_dirDataGCamp(dLGC,1).name;
    
    for an = AnList
        if strncmpi(an{1},NameAnimal,length(an{1})) > 0
            
            NameAnimal
            out_cell{row,1}=NameAnimal;
            PathFile_GCamp = [dirAllDataGCamp,filesep,NameAnimal];
            dirNameGCamp   = dir(PathFile_GCamp);
            
            for dNGC=3:length(dirNameGCamp)
                NameFile_GCamp_Day =  dirNameGCamp(dNGC,1).name;
                if contains(NameFile_GCamp_Day,'_Par_MIPSIP_Par') %if exist
                    
                    %load GCamp
                    load([PathFile_GCamp,filesep,NameFile_GCamp_Day])
                    DayFile = dataGCamp.Info.Date;
                    out_cell{row,str2num(dataGCamp.Info.Date)+1}=dataGCamp.InfoTrial.NumTrials;
                    
                end
            end
            
            row=row+1;
        end
    end
end