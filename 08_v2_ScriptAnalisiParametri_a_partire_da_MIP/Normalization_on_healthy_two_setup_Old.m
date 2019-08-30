% NORMALIZATION on healthy and sani for different setup

%        InfoPicchi: column 1 ->
%                   column 2 ->
%                   column 3 -> peak start point
%                   column 4 -> peak duration
%                   column 5 -> max value peak
%                   column 6 -> max value peak point (IstAmpMaxWSig)
%                   column 7 -> full width at half maximum (FWHM)
%                   column 8 -> area under the peak curve (AUC)
%                   column 9 -> peak-to-peak amplitude (PtPAmp)
%                   column 10 -> NumberOfRealPeak inside the peak     %non più il num of subpeaks inside the main peak (numPks)
%                   column 11 -> SlopeInitial (slope between onset _ max value)
%                   column 12 -> Time to Peak ( [max value peak point - peak start point] )
%                   column 13 -> Time Distance bw first movement - onset
%                   column 14 -> Time Distance bw first movement - time of max


%last updating 27 04 2016

%add different values of treats in the new setup if other groups
%%% number ID for each group
% control = 0; stroke = 1; robot+Bont_4w = 2; robot+Bont_1w = 3;
% robot_4w = 21; robot_1w = 31; robot+stim_4w=214; TOX =5; 6 sani_new
% setup; 7 robot_4w_new setup
clear all

old_setup=[0 1 2 3 21 31];
new_setup=[5 6 7 8 9 10 11 12]; %sostitution of the number of treatment (6->0 and 7->21, to use the normal code)

name='dataMouseGCamp_MEAN_SelROI_1_02-May-2018_12-9';
name='dataMouseGCamp_MEAN_SelROI_1_15-Jun-2018_10-33';
load(name)

if ~exist('MeanPar')
    error('LOAD the matrix "dataMouseGCamp_Median_Sel....mat" ')
end

%% output data for statistical analysis

% Ale 180420
% assuming that the data are the normalized one...



par_list = [11 8 5 6];
par_name = {'Slope', 'AUC', 'Peak', 'DeltaT'};
group_list = [0 1 3 2 21 ...
    6 7 5 8 ...
    12 10 11 9 ...
    31];
group_name = {'ctrl','strk', 'strk-bont', 'reh4w', 'robot',...
    'ctrl_{new}','robot_{new}','tox_{new}','toxA_{new}', ...
    'sham',  'optostim', 'robot_{or}', 'optostim+robot', ...
    'robot'};

animal_list = {};
%creating list of animals based on groups
i_par = 11;
for i_day=1:size(MeanPar,2)
    MeanPar_par_one=MeanPar{i_par,i_day};
    for i_mice=1:size(MeanPar_par_one,1)
        animal_name = MeanPar_par_one{i_mice,2};
        tr=MeanPar_par_one{i_mice,4};
        if ismember (tr,group_list)
            if ~isempty(animal_name)
                if ~ismember(animal_name, animal_list)
                    animal_list = cat(2,animal_list, {animal_name});
                end
            end
        end
    end
end

output = {};
header = {};
for i_par=par_list
    disp(['Analyzing parameter:',num2str(i_par)])
    
    for i_day=1:size(MeanPar,2)
        disp(['Extracting day:',num2str(i_day)])
        col_num = (find(par_list==i_par)-1)*size(MeanPar,2) + i_day;
        MeanPar_par_one=MeanPar{i_par,i_day};

        header{col_num}=[par_name{find(par_list==i_par)},'_D',num2str(i_day)];
        header{length(par_list)*size(MeanPar,2) + 1} = 'group_name';
        header{length(par_list)*size(MeanPar,2) + 2} = 'subject';

        for i_mice=1:size(MeanPar_par_one,1)
            
            
            tr=MeanPar_par_one{i_mice,4};
            if ismember (tr,group_list)
                row = find(strcmpi(MeanPar_par_one{i_mice,2}, animal_list));
                
                disp(['Extracting subject:',MeanPar_par_one{i_mice,2}])
                output{row, col_num} = MeanPar_par_one{i_mice,1}(2,3);
                output{row, length(par_list)*size(MeanPar,2) + 1} = group_name{find(group_list==tr)};
                output{row, length(par_list)*size(MeanPar,2) + 2} = MeanPar_par_one{i_mice,2};
            end
        end
        MeanPar{i_par,i_day}=MeanPar_par_one;
    end
end
output= cat(1, header, output);
tableout = cell2table(output(2:end,:),'VariableNames',output(1,:));
writetable(tableout,[name,'.csv'])

%% Find mean parameters of healthy and sani groups
Mean_healthy=zeros(14,3,2);
Mean_sani=zeros(14,3,2);

for i_par=3:size(MeanPar,1)
    all_data_healthy_fluo=[];
    all_data_healthy_force=[];
    all_data_sani_fluo=[];
    all_data_sani_force=[];
    for i_day=1:size(MeanPar,2)
        MeanPar_par_one=MeanPar{i_par,i_day};
        for i_mice=1:size(MeanPar_par_one,1)
            if MeanPar_par_one{i_mice,4}==0
                all_data_healthy_fluo=[all_data_healthy_fluo, MeanPar_par_one{i_mice,1}(:,3)];
                all_data_healthy_force=[all_data_healthy_force, MeanPar_par_one{i_mice,1}(:,1)];
            elseif MeanPar_par_one{i_mice,4}==6
                all_data_sani_fluo=[all_data_sani_fluo, MeanPar_par_one{i_mice,1}(:,3)];
                all_data_sani_force=[all_data_sani_force, MeanPar_par_one{i_mice,1}(:,1)];
            end
        end
    end
    Mean_healthy(i_par,:,1)=nanmean(all_data_healthy_force,2);
    Mean_healthy(i_par,:,2)=nanmean(all_data_healthy_fluo,2);
    Mean_sani(i_par,:,1)=nanmean(all_data_sani_force,2);
    Mean_sani(i_par,:,2)=nanmean(all_data_sani_fluo,2);
end

clearvars -except Mean_healthy Mean_sani name MeanPar old_setup new_setup

%% Normalization of all parameters

for i_par=3:size(MeanPar,1)
    for i_day=1:size(MeanPar,2)
        MeanPar_par_one=MeanPar{i_par,i_day};
        for i_mice=1:size(MeanPar_par_one,1)
            tr=MeanPar_par_one{i_mice,4};
            if ismember (tr,old_setup)
                MeanPar_par_one{i_mice,1}(:,3)=MeanPar_par_one{i_mice,1}(:,3)./Mean_healthy(i_par,:,2)';
                MeanPar_par_one{i_mice,1}(:,1)=MeanPar_par_one{i_mice,1}(:,1)./Mean_healthy(i_par,:,1)';
            elseif ismember(tr,new_setup)
                MeanPar_par_one{i_mice,1}(:,3)=MeanPar_par_one{i_mice,1}(:,3)./Mean_sani(i_par,:,2)';
                MeanPar_par_one{i_mice,1}(:,1)=MeanPar_par_one{i_mice,1}(:,1)./Mean_sani(i_par,:,1)';
                %                 if tr==6
                %                     MeanPar_par_one{i_mice,4}=0;
                %                 elseif tr==7
                %                     MeanPar_par_one{i_mice,4}=21;
                %                 end
            elseif isempty(tr)
                disp('missing data')
            else
                disp(tr)
                disp('group not identified')
            end
        end
        MeanPar{i_par,i_day}=MeanPar_par_one;
    end
end

Filename=[name(1:end-5),'_normalized'];
save(Filename,'MeanPar')


%% output data for statistical analysis

% Ale 180420
% assuming that the data are the normalized one...



par_list = [11 8 5 6];
par_name = {'Slope', 'AUC', 'Peak', 'DeltaT'};
group_list = [0 1 3 2 21 ...
    6 7 5 8 ...
    12 10 11 9 ...
    31];
group_name = {'ctrl','strk', 'strk-bont', 'reh4w', 'robot',...
    'ctrl_{new}','robot_{new}','tox_{new}','toxA_{new}', ...
    'sham',  'optostim', 'robot_{or}', 'optostim+robot', ...
    'robot'};

animal_list = {};
%creating list of animals based on groups
i_par = 11;
for i_day=1:size(MeanPar,2)
    MeanPar_par_one=MeanPar{i_par,i_day};
    for i_mice=1:size(MeanPar_par_one,1)
        animal_name = MeanPar_par_one{i_mice,2};
        tr=MeanPar_par_one{i_mice,4};
        if ismember (tr,group_list)
            if ~isempty(animal_name)
                if ~ismember(animal_name, animal_list)
                    animal_list = cat(2,animal_list, {animal_name});
                end
            end
        end
    end
end

output = {};
header = {};
for i_par=par_list
    disp(['Analyzing parameter:',num2str(i_par)])
    
    for i_day=1:size(MeanPar,2)
        disp(['Extracting day:',num2str(i_day)])
        col_num = (find(par_list==i_par)-1)*size(MeanPar,2) + i_day;
        MeanPar_par_one=MeanPar{i_par,i_day};
        
        header{col_num}=[par_name{find(par_list==i_par)},'_D',num2str(i_day)];
        header{length(par_list)*size(MeanPar,2) + 1} = 'group_name';
        header{length(par_list)*size(MeanPar,2) + 2} = 'subject';
        
        
        for i_mice=1:size(MeanPar_par_one,1)
            
            
            tr=MeanPar_par_one{i_mice,4};
            if ismember (tr,group_list)
                row = find(strcmpi(MeanPar_par_one{i_mice,2}, animal_list));
                
                disp(['Extracting subject:',MeanPar_par_one{i_mice,2}])
                output{row, col_num} = MeanPar_par_one{i_mice,1}(2,3);
                output{row, length(par_list)*size(MeanPar,2) + 1} = group_name{find(group_list==tr)};
                output{row, length(par_list)*size(MeanPar,2) + 2} = MeanPar_par_one{i_mice,2};
            end
        end
        MeanPar{i_par,i_day}=MeanPar_par_one;
    end
end
output= cat(1, header, output);
tableout = cell2table(output(2:end,:),'VariableNames',output(1,:));
writetable(tableout,[name,'_normalized.csv'])
