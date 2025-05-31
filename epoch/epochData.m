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
%
% example inputs:
% subject = 'NU002';
% indir = '/project/joelvoss/data/';
% outdir = '/project/joelvoss/tmp-rostowsky/testOuts';
% codedir = '/project/joelvoss/code/';
% pathToImages = '/project/joelvoss/code/retinotopy/workspace_retinotopyCaltsmash.mat';
%
% Output
% epochData.mat
% outdir/logs/epochData_log.txt

% set up paths
addpath(genpath('/project/joelvoss/tmp-rostowsky/github/retinotopy/dependencies/epoch')); % this has the get_dirs_kr function
addpath(genpath([codedir '/retinotopy' ]));
addpath(genpath([codedir '/Nlx2Mat' ]));

% load behavioral data for subject
behavFiles = load_behavioral_data(indir, subject);
retinotopyImages = load(pathToImages);

% find data and image directory for subject, load MNI
[csc_dir, img_dir] = get_dirs_kr(indir, subject);

% loads depth electrode labels
findDepthElectrodes = readtable([img_dir '/' subject '_mni.csv']).channel;

% unique depth electrodes
[uniqueDepth] = unique_depth_electrodes(findDepthElectrodes);

% organized depth electrodes
contacts = readtable([img_dir 'jackbox.csv']).Elect; % the jacksheet contains the correct order of contacts
findDepthElectrodes = reorder_depth_electrodes(findDepthElectrodes, uniqueDepth, contacts);

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
save([outdir '/' 'epochData.mat']);
fprintf('\nfinished\n');

end
