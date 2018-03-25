%%%%%%% Info Animal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [num_week num_day_per_week list_days] = Fun_Info_Animal(Animal_name,list_real_days);
function [weekly_day_stop_All , num_week] = Fun_Info_Animal(Animal_name,list_real_days)

if strfind(Animal_name,'control')
    %control
    num_week         = 1;
    num_day_per_week = 5;
    list_days = [1:5];
    
elseif strfind(Animal_name,'stroke')
    
    if strfind(Animal_name,'BoNT')
        %rehab
        num_week         = 4;
        num_day_per_week = 5;
        %list real day
        list_days = [1:5;6:10;11:15;16:20];
    else
        %stroke
        num_week         = 1;
        num_day_per_week = 5;
        list_days = [1:5];
        %exception
        if strcmp(Animal_name,'GCaMPChR2_25_stroke')
            num_day_per_week = 4;
            list_days = 1:4;
        end
    end
end

%find last day of every week
weekly_day_stop_All = [];
for i_w=1:num_week
    
    for i_d=1:  num_day_per_week
        weekly_day_stop = [find(list_real_days == list_days(i_w,num_day_per_week-i_d+1))];
        if ~isempty(weekly_day_stop)
            weekly_day_stop_All = [weekly_day_stop_All; list_real_days(weekly_day_stop)];
            break
        end
    end   
end
