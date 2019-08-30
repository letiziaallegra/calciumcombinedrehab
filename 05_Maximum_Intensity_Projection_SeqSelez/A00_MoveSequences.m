function status = A00_MoveSequences(main_dir)

% making sure it ends with a fileseparator
if main_dir(end) ~= filesep
    main_dir = [main_dir,filesep];
end

% this is where the IMAGE SEQUENCE ROOT DIR SHOULD BE
WORKING_FOLDER = [fileparts(main_dir(1:end-1)),filesep,'Sequence_Trials_DATA/'];

if ~isdir(WORKING_FOLDER)
    error([WORKING_FOLDER, ' does not exist'])
end

if WORKING_FOLDER(end) ~= filesep
    WORKING_FOLDER = [main_dir,filesep];
end

sub_folders = dir(main_dir);
for subidx = 3:length(sub_folders)
    sub_folder = sub_folders(subidx);
    if ~ sub_folder.isdir
        continue
    end
    day_folders = dir([sub_folder.folder,filesep,sub_folder.name]);
    for dayidx = 3:length(day_folders)
        day_folder = day_folders(dayidx);
        if ~ day_folder.isdir
            continue
        end
        seq_folder_name = [sub_folder.name,'_',day_folder.name(1:2),'_SequenceLong_TIF',filesep,'SequenceLong',filesep];
        src_folder = [WORKING_FOLDER,sub_folder.name,filesep,seq_folder_name];
        dst_folder = [day_folder.folder,filesep,day_folder.name,filesep,'SequenceLong',filesep];
        try
            copyfile(src_folder,dst_folder,'f')
        catch
            disp([dst_folder, ' not a valid day or sub folder'])
        end
    end
    
end