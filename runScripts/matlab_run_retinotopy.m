clc; clear;
subjectDir = '/project/joelvoss/tmp-rostowsky/hpcData';
indir = '/project/joelvoss/data';
outdir = '/project/joelvoss/tmp-rostowsky/hpcData';
codedir = '/project/joelvoss/code';

pathToImages = '/project/joelvoss/code/retinotopy/workspace_retinotopyCaltsmash.mat';
subjects = struct2cell(dir(subjectDir));
subjects = subjects(1,3:end);

%% step 1
for j = 1:length(subjects)
    epochData(subjects{j}, indir,[outdir '/' subjects{j}], codedir, pathToImages);
end

%%
figdir = '/project/joelvoss/tmp-rostowsky/hpcDataFigs';
dlDir = '/project/joelvoss/tmp-rostowsky/dlDir'; % this is directory that you can download onto your local machine to look at the images 
for j = 1:length(subjects)
    mkdir([dlDir '/' subjects{j} '/epoched/' ]);
    [resultLabel] = matchSubjectToHpcElectrodes(subjects{j});
    currSubjectFigures = struct2cell(dir([figdir '/' subjects{j}]));
    currSubjectFigures = currSubjectFigures(1,3:end);
    for k = 1:length(resultLabel)
        fileName = find(contains(currSubjectFigures, ['-' resultLabel{k} '.png']));
        if length(fileName) > 1
            error('too many files matched, try regexp instead');
        end
        copyfile([figdir '/' subjects{j} '/' currSubjectFigures{fileName}], [dlDir '/' subjects{j} '/epoched/']);
    end
end
for j = 1:length(subjects)
    mkdir([dlDir '/' subjects{j} '/behavioral/' ]);
    copyfile([subjectDir '/' subjects{j} '/*.png'], [dlDir '/' subjects{j} '/behavioral/' ]);
end

%% step 2
parfor j = 1:length(subjects)
    rereferenceElectrodes([subjectDir '/' subjects{j}]);
end

%% detect spikes, don't remove yet
parfor j = 1:length(subjects)
    removeEpileptiformActivity([subjectDir '/' subjects{j}]);
end

%% select data from hippocampal electrodes & essentially resample stimulus to the sampling frequency
parfor j = 1:length(subjects)
    extractHpcData([subjectDir '/' subjects{j}], subjects{j})
end

%% run TFA
parfor j = 1:length(subjects)
    tfAnalysis_jk_kr(subjectDir, subjects{j});
end

%%