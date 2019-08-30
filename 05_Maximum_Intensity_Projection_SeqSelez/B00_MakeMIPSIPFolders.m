function status = B00_MakeMIPSIPFolders(main_dir, store_folder)

% making sure it ends with a fileseparator
if main_dir(end) ~= filesep
    main_dir = [main_dir,filesep];
end

% making sure it ends with a fileseparator
if store_folder(end) ~= filesep
    store_folder = [store_folder,filesep];
end


src_folder = main_dir;
[~, animal_name] = fileparts(main_dir(1:end-1));
dst_folder = [store_folder,animal_name];
copyfile([src_folder,'*.mat'],dst_folder,'f')
copyfile([src_folder,'*.fig'],dst_folder,'f')
copyfile([src_folder,'*.tif'],dst_folder,'f')

