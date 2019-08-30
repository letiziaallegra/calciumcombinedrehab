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

old_setup=[0 1 2 3 21 31 100 101];
new_setup=[5 6 7 8 9 10 11 12 61]; %sostitution of the number of treatment (6->0 and 7->21, to use the normal code)

name='dataMouseGCamp_MEAN_SelROI_1_02-May-2018_12-9';
name='dataMouseGCamp_MEAN_SelROI_1_15-Jun-2018_10-33';
name='dataMouseGCamp_MEAN_SelROI_1_22-Jun-2018_10-53';
name='dataMouseGCamp_MEAN_SelROI_1_08-Oct-2018_16-39';
name='dataMouseGCamp_MEAN_SelROI_1_22-Jan-2019_14-46';
name='dataMouseGCamp_MEAN_SelROI_1_23-Jan-2019_11-7';
name='dataMouseGCamp_MEAN_SelROI_1_23-Jan-2019_13-34';

load(name)

if ~exist('MeanPar')
    error('LOAD the matrix "dataMouseGCamp_Median_Sel....mat" ')
end

%% output data for statistical analysis

% Ale 180420
% assuming that the data are the normalized one...



par_list = [11 8 5 6 61 66];
par_name = {'Slope', 'AUC', 'Peak', 'TimeToPeakRostral', ...
    'TimeToPeakCaudal', 'DeltaT'};
group_list = [0 100 101 1 3 2 21 ...
    6 7 5 8 ...
    12 10 11 9 ...
    31 61];
group_name = {'ctrl', 'new-ctrl-paper', 'new-ctrl-paper-1', 'strk', 'strk-bont', 'reh4w', 'robot',...
    'ctrl_{new}','robot_{new}','tox_{new}','toxA_{new}', ...
    'sham',  'optostim', 'robot_{or}', 'optostim+robot', ...
    'robot','ctrl_{new}^{long}'};

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
    par_start = (find(par_list==i_par)-1)*size(MeanPar,2);
    header_name = par_name{find(par_list==i_par)};
    disp(['Analyzing parameter:',num2str(i_par)])
    if i_par == 61
        i_par = 6;
        prow = 3;
    elseif i_par == 66
        i_par = 6;
        prow = [3,2];
    else
        prow = 2;
    end
    for i_day=1:size(MeanPar,2)

        disp(['Extracting day:',num2str(i_day)])
        col_num = par_start + i_day;
        MeanPar_par_one=MeanPar{i_par,i_day};

        header{col_num}=[header_name,'_D',num2str(i_day)];
        header{length(par_list)*size(MeanPar,2) + 1} = 'group_name';
        header{length(par_list)*size(MeanPar,2) + 2} = 'subject';
        header{length(par_list)*size(MeanPar,2) + 3} = 'group_number';

        for i_mice=1:size(MeanPar_par_one,1)
            
            
            tr=MeanPar_par_one{i_mice,4};
            if ismember (tr,group_list)
                row = find(strcmpi(MeanPar_par_one{i_mice,2}, animal_list));
                disp(['Extracting subject:',MeanPar_par_one{i_mice,2}])
                if length(prow)==1
                    output{row, col_num} = MeanPar_par_one{i_mice,1}(prow,3);
                else
                    output{row, col_num} = MeanPar_par_one{i_mice,1}(prow(1),3) - MeanPar_par_one{i_mice,1}(prow(2),3);
                end
                output{row, length(par_list)*size(MeanPar,2) + 1} = group_name{find(group_list==tr)};
                output{row, length(par_list)*size(MeanPar,2) + 2} = MeanPar_par_one{i_mice,2};
                output{row, length(par_list)*size(MeanPar,2) + 3} = tr;
            end
        end
        MeanPar{i_par,i_day}=MeanPar_par_one;
    end
end
output= cat(1, header, output);
tableout = cell2table(output(2:end,:),'VariableNames',output(1,:));
writetable(tableout,[name,'.csv'])
%% normalization

for col = 1:size(tableout,2) - 3
    % normalizza per giorno
    % mean_old = nanmean(cell2mat(tableout{tableout{:,33}==0,col}));
    % mean_new = nanmean(cell2mat(tableout{tableout{:,33}==6,col}));
    % normalizza per settimana
    if mod(col,5) == 1
        mean_old = nanmean(cell2mat(reshape(tableout{tableout{:,33}==0,(ceil(col/5)-1)*5+1:ceil(col/5)*5},1,[])));
        mean_new = nanmean(cell2mat(reshape(tableout{tableout{:,33}==6,(ceil(col/5)-1)*5+1:ceil(col/5)*5},1,[])));
    end
    
    for row = 1:size(tableout,1)
        if ismember(tableout{row,33}, old_setup)
            tableout{row,col} = {cell2mat(tableout{row,col}) / mean_old};
        else
            tableout{row,col} = {cell2mat(tableout{row,col}) / mean_new};
        end
    end
end
writetable(tableout,[name,'_normalized.csv'])
