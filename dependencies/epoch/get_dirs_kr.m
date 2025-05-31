function [csc_dir, img_dir] = get_dirs_kr(indir,subject)
% gets the directory information for a given subject

% this function really is just a way of keeping track of where the data for
% this experiment lives in session directories


checkDirs = struct2cell(dir([indir '/' subject '/ieeg' ]));
checkDirs = checkDirs(1,3:end);
regex = '[0-9]{4}-[0-9]{2}-[0-9]{2}';
dateFind = regexp(checkDirs, regex);
for j = 1:length(dateFind)
    if isempty(dateFind{j})
        dateFind{j} = 0;
    end
end
dateFind = checkDirs(find([dateFind{:}] ~= 0));

for j = length(dateFind):-1:1
    if size(dir([indir '/' subject '/ieeg/' checkDirs{j}]),1) == 2
        dateFind(j) = [];
    end
end

if length(dateFind) ~= 1 && ~isempty(dateFind)
    compDate = struct2cell(dir([indir '/' subject '/behav/retinotopy']));
    compDate = compDate(1,3:end);
    compDate = strsplit(compDate{1},'_');
    compDate = compDate{1}(1:8);
    compDate = datetime(compDate,'InputFormat', 'yyyyMMdd');

    for j = length(dateFind):-1:1
        currDate = strsplit(dateFind{j},'_');
        currDate = datetime(currDate{1}, 'InputFormat', 'yyyy-MM-dd');
        if currDate ~= compDate
            dateFind(j) = [];
        end
    end
end

if isempty(dateFind)
    checkEdf = dir([indir '/' subject '/ieeg/*RETINOTOPY.edf']);
else
    checkEdf = [];
end

if length(dateFind) ~= 1 && isempty(checkEdf)
    error('Cannot continue, no matching iEEG and behavior dates found');
end

if ~isempty(checkEdf)
    csc_dir = ['/project/joelvoss/data/' subject '/ieeg/']; 
else
    sess_dir = dateFind{1};
    csc_dir = ['/project/joelvoss/data/' subject '/ieeg/' sess_dir filesep]; %TODO: which session corresponds to a BORG datafile
end
img_dir = ['/project/joelvoss/data/' subject '/anat/'];
end


