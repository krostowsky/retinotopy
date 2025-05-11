function epochData(subject, indir, outdir, codedir, pathToImages)

% Description: epochData aligns data from separate behavioral trials (blocks) used
% in retinotopy experiments
%
% epochData(subject, indir, outdir, codedir, pathToImages)
% 
% Input
% - subject   = 1xN character array with subject id (e.g. UC006) in indir
% - indir     = directory where subject data is held
% - outdir    = directory where data will be saved to
% - codedir   = directory with neuralynx utilities and epoch scripts
% - pathToImages  = full path to .mat file with images used in experiment

% example inputs:
% subject = 'NU002';
% indir = '/project/joelvoss/data/';
% outdir = '/project/joelvoss/tmp-rostowsky/testOuts';
% codedir = '/project/joelvoss/code/';
% pathToImages = '/project/joelvoss/code/retinotopy/workspace_retinotopyCaltsmash.mat';


%% set up paths
addpath(genpath('/scratch/midway3/rostowsky/gitlab/retinotopy_kr_code/code/projectSpecificCode/epochScripts')); % this has the get_dirs_kr function
addpath(genpath([codedir '/retinotopy' ]));
addpath(genpath([codedir '/Nlx2Mat' ]));

%% load behavioral data for subject
behavFiles = struct2cell(dir([indir '/' subject '/behav/retinotopy/*.mat']));
behavFiles = behavFiles(1,:);
retinotopyImages = load(pathToImages);

%% find data and image directory for subject, load MNI
[csc_dir, img_dir] = get_dirs_kr(indir, subject);

% the mni file contains ONLY depth electrode data, this is a quick way
% to reject scalp electrodes
findDepthElectrodes = readtable([img_dir '/' subject '_mni.csv']);
findDepthElectrodes = findDepthElectrodes.channel;

%% organize depth electrodes from MNI
organizeDepthElectrodes = cell(length(findDepthElectrodes),1);
for j = 1:length(findDepthElectrodes)
    charInd = regexp(findDepthElectrodes{j}, '[A-z]');
    organizeDepthElectrodes{j} = findDepthElectrodes{j}(charInd(1):charInd(end));
end

uniqueDepth = unique(organizeDepthElectrodes,'stable');
countElectrode = 0;
for j = 1:length(uniqueDepth)
    currMatch = find(strcmpi(uniqueDepth{j}, organizeDepthElectrodes));
    for k = 1:length(currMatch)
        countElectrode = countElectrode + 1;
        organizeDepthElectrodes{countElectrode} = [uniqueDepth{j} num2str(k)];
    end
end

jacksheet = readtable([img_dir 'jackbox.csv']);
contacts = jacksheet.Elect;
orderElect = cell(length(uniqueDepth),2); % first column is the depth electrode, second column is the jackbox (CSC file?) location
contactNew = cell(length(contacts),1);
contactNames = regexp(contacts, '[A-z]');
for j = 1:length(contactNew)
    contactNew{j} = contacts{j}(contactNames{j});
end

for j = 1:length(orderElect)
    orderElect{j,1} = uniqueDepth{j};
    currFind = find(strcmpi(uniqueDepth{j},contactNew));
    orderElect{j,2} = currFind(1); % first index with electrode match
end

%% this orders the depth electrodes by their appearance in the jacksheet 
tmpElects = cell(size(findDepthElectrodes,1),1);
[~, sortOrder] = sort(cell2mat(orderElect(:,2)));
orderElect(:,2) = num2cell(sortOrder); % second column is the row index of the replacement row
currCountElect = 0;
for j = 1:size(orderElect,1)
    currLabel = orderElect{orderElect{j,2},1};
    currFind = find(strcmpi(currLabel, contactNew));
    currResult = arrayfun(@(x) sprintf('%s%d', currLabel, x), 1:length(currFind), 'UniformOutput', false); % assign number for each of the current label
    tmpElects(currCountElect+1:currCountElect+length(currFind)) = currResult;
    currCountElect = currCountElect + length(currFind);
end

findDepthElectrodes = tmpElects;
clearvars -except subject csc_dir img_dir behavFiles indir codedir retinotopyImages findDepthElectrodes outdir

%%
[datCell, outStruct, rawOnOffInds, rawOnOffId, sortedOnOffInds, sortedOnOffId, builtTimes, corruptElectrodes, sortedBlockInds, sortedRawOnOffId, sortedRawOnOffInds] = get_subj_eeg_kr(subject, ...
    csc_dir, ...
    img_dir, ...
    behavFiles, ...
    indir, ...
    codedir, ...
    retinotopyImages, ...
    findDepthElectrodes, outdir);

%%
%findDepthElectrodes(logical(corruptElectrodes)) = [];

%%
save([outdir '/' 'epochData.mat']);
fprintf('\nfinished\n');

end
