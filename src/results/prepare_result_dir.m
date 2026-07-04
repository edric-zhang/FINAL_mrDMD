function prepare_result_dir(result_dir, clear_existing)
%PREPARE_RESULT_DIR Create a result directory and optionally clear old files.

if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

if clear_existing
    files = dir(result_dir);
    for k = 1:numel(files)
        if files(k).isdir
            continue;
        end
        delete(fullfile(files(k).folder, files(k).name));
    end
end

end
